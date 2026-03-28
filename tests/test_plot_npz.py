"""
tests/test_plot_npz.py – unit tests for pyaxon/plot_npz.py.

Each TestCase class covers one public function.  Matplotlib is switched to the
non-interactive 'Agg' backend so no windows ever open during the test run.
All file I/O uses ``tempfile.TemporaryDirectory`` so every test is isolated.
"""

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")   # must be set before any other matplotlib import
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

from pyaxon.plot_npz import (
    DEFAULT_FIELDS,
    _build_parser,
    _step_sort_key,
    collect_npz_files,
    load_simulation_npz,
    make_gif_from_npz_files,
    plot_npz,
    plot_simulation_snapshot,
)
from tests.fixtures import (
    SHAPE,
    make_snapshot,
    make_snapshot_dir,
    make_snapshot_no_vm,
    make_snapshot_partial,
)


# =============================================================================
# load_simulation_npz
# =============================================================================
class TestLoadSimulationNpz(unittest.TestCase):

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)

    def tearDown(self):
        plt.close("all")
        self._tmp.cleanup()

    def test_returns_dict_with_all_fields(self):
        p = make_snapshot(self.tmp / "snap.npz")
        data = load_simulation_npz(p)
        for field in DEFAULT_FIELDS:
            self.assertIn(field, data)
        self.assertIn("v_m", data)

    def test_values_are_ndarrays(self):
        p = make_snapshot(self.tmp / "snap.npz")
        data = load_simulation_npz(p)
        for v in data.values():
            self.assertIsInstance(v, np.ndarray)

    def test_scalar_field_shape_matches(self):
        p = make_snapshot(self.tmp / "snap.npz", shape=SHAPE)
        data = load_simulation_npz(p)
        self.assertEqual(data["phi"].shape, SHAPE)

    def test_vm_field_shape_is_2_x_nx_x_ny(self):
        p = make_snapshot(self.tmp / "snap.npz", shape=SHAPE)
        data = load_simulation_npz(p)
        self.assertEqual(data["v_m"].shape, (2,) + SHAPE)

    def test_raises_file_not_found_for_missing_file(self):
        with self.assertRaises(FileNotFoundError):
            load_simulation_npz(self.tmp / "nonexistent.npz")

    def test_loads_partial_snapshot(self):
        p = make_snapshot_partial(self.tmp / "partial.npz", ["phi", "psi"])
        data = load_simulation_npz(p)
        self.assertIn("phi", data)
        self.assertIn("psi", data)
        self.assertNotIn("ngf", data)

    def test_loads_snapshot_without_vm(self):
        p = make_snapshot_no_vm(self.tmp / "no_vm.npz")
        data = load_simulation_npz(p)
        self.assertNotIn("v_m", data)

    def test_string_path_accepted(self):
        p = make_snapshot(self.tmp / "snap.npz")
        data = load_simulation_npz(str(p))
        self.assertIn("phi", data)


# =============================================================================
# _step_sort_key
# =============================================================================
class TestStepSortKey(unittest.TestCase):

    def test_trailing_number_is_used(self):
        key = _step_sort_key(Path("sim_10000.npz"))
        self.assertEqual(key[0], 10000)

    def test_larger_step_sorts_later(self):
        k1 = _step_sort_key(Path("sim_0.npz"))
        k2 = _step_sort_key(Path("sim_10000.npz"))
        k3 = _step_sort_key(Path("sim_200000.npz"))
        self.assertLess(k1, k2)
        self.assertLess(k2, k3)

    def test_no_number_returns_inf(self):
        key = _step_sort_key(Path("snapshot.npz"))
        self.assertEqual(key[0], float("inf"))

    def test_files_sorted_chronologically(self):
        names = ["sim_20000.npz", "sim_0.npz", "sim_10000.npz"]
        paths = [Path(n) for n in names]
        sorted_paths = sorted(paths, key=_step_sort_key)
        self.assertEqual([p.name for p in sorted_paths],
                         ["sim_0.npz", "sim_10000.npz", "sim_20000.npz"])


# =============================================================================
# collect_npz_files
# =============================================================================
class TestCollectNpzFiles(unittest.TestCase):

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def test_single_file_returns_list_of_one(self):
        p = make_snapshot(self.tmp / "snap.npz")
        result = collect_npz_files(p)
        self.assertEqual(result, [p])

    def test_directory_returns_sorted_files(self):
        paths = make_snapshot_dir(self.tmp / "run", n_steps=3, prefix="s_")
        result = collect_npz_files(self.tmp / "run")
        self.assertEqual(result, paths)

    def test_directory_order_is_chronological(self):
        make_snapshot_dir(self.tmp / "run", n_steps=4, step_size=10000, prefix="s_")
        result = collect_npz_files(self.tmp / "run")
        steps = [int(p.stem.split("_")[1]) for p in result]
        self.assertEqual(steps, sorted(steps))

    def test_raises_for_empty_directory(self):
        empty = self.tmp / "empty"
        empty.mkdir()
        with self.assertRaises(ValueError):
            collect_npz_files(empty)

    def test_glob_pattern_returns_matching_files(self):
        make_snapshot_dir(self.tmp / "run", n_steps=3, prefix="g_")
        result = collect_npz_files(self.tmp / "run" / "g_*.npz")
        self.assertEqual(len(result), 3)

    def test_raises_for_no_glob_matches(self):
        with self.assertRaises(ValueError):
            collect_npz_files(self.tmp / "run" / "*.npz")

    def test_string_path_accepted(self):
        p = make_snapshot(self.tmp / "snap.npz")
        result = collect_npz_files(str(p))
        self.assertEqual(len(result), 1)


# =============================================================================
# plot_simulation_snapshot
# =============================================================================
class TestPlotSimulationSnapshot(unittest.TestCase):

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)

    def tearDown(self):
        plt.close("all")
        self._tmp.cleanup()

    def _load(self, include_vm=True):
        p = make_snapshot(self.tmp / "snap.npz", include_vm=include_vm)
        return load_simulation_npz(p)

    def test_returns_figure_and_axes(self):
        data = self._load()
        fig, axes = plot_simulation_snapshot(data)
        self.assertIsInstance(fig, plt.Figure)
        self.assertEqual(len(axes), 8)

    def test_title_is_set_when_provided(self):
        data = self._load()
        fig, _ = plot_simulation_snapshot(data, title="test title")
        self.assertEqual(fig.texts[0].get_text(), "test title")

    def test_no_title_when_not_provided(self):
        data = self._load()
        fig, _ = plot_simulation_snapshot(data, title=None)
        self.assertEqual(len(fig.texts), 0)

    def test_handles_missing_vm(self):
        data = self._load(include_vm=False)
        fig, axes = plot_simulation_snapshot(data)
        self.assertIsInstance(fig, plt.Figure)

    def test_handles_partial_fields(self):
        p = make_snapshot_partial(self.tmp / "p.npz", ["phi", "psi"])
        data = load_simulation_npz(p)
        fig, axes = plot_simulation_snapshot(data)
        self.assertIsInstance(fig, plt.Figure)

    def test_cmap_parameter_is_accepted(self):
        data = self._load()
        fig, _ = plot_simulation_snapshot(data, cmap="plasma")
        self.assertIsInstance(fig, plt.Figure)

    def test_vm_2d_fallback_path(self):
        """v_m with unexpected shape should still render without error."""
        p = make_snapshot(self.tmp / "snap.npz",
                          extra_fields={"v_m": np.ones(SHAPE)})
        data = load_simulation_npz(p)
        fig, _ = plot_simulation_snapshot(data)
        self.assertIsInstance(fig, plt.Figure)

    def test_unused_axes_are_turned_off(self):
        """If fewer than 8 fields are present, remaining axes should be invisible."""
        p = make_snapshot_partial(self.tmp / "p.npz", ["phi"])
        data = load_simulation_npz(p)
        fig, axes = plot_simulation_snapshot(data)
        invisible = [not ax.get_visible() or not ax.axison for ax in axes]
        # At least some axes must be invisible / disabled.
        self.assertTrue(any(invisible))


# =============================================================================
# plot_npz
# =============================================================================
class TestPlotNpz(unittest.TestCase):

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)

    def tearDown(self):
        plt.close("all")
        self._tmp.cleanup()

    def test_returns_figure_and_axes(self):
        p = make_snapshot(self.tmp / "snap.npz")
        fig, axes = plot_npz(p, show=False)
        self.assertIsInstance(fig, plt.Figure)

    def test_saves_image_when_save_path_provided(self):
        p = make_snapshot(self.tmp / "snap.npz")
        out = self.tmp / "out.png"
        plot_npz(p, save_path=out, show=False)
        self.assertTrue(out.exists())
        self.assertGreater(out.stat().st_size, 0)

    def test_does_not_save_when_save_path_is_none(self):
        p = make_snapshot(self.tmp / "snap.npz")
        plot_npz(p, save_path=None, show=False)
        png_files = list(self.tmp.glob("*.png"))
        self.assertEqual(png_files, [])

    def test_show_false_does_not_call_plt_show(self):
        p = make_snapshot(self.tmp / "snap.npz")
        with patch("pyaxon.plot_npz.plt.show") as mock_show:
            plot_npz(p, show=False)
            mock_show.assert_not_called()

    def test_show_true_calls_plt_show(self):
        p = make_snapshot(self.tmp / "snap.npz")
        with patch("pyaxon.plot_npz.plt.show") as mock_show:
            plot_npz(p, show=True)
            mock_show.assert_called_once()

    def test_raises_for_missing_file(self):
        with self.assertRaises(FileNotFoundError):
            plot_npz(self.tmp / "missing.npz", show=False)

    def test_cmap_forwarded_to_snapshot(self):
        p = make_snapshot(self.tmp / "snap.npz")
        # Should not raise; cmap is forwarded without validation.
        fig, _ = plot_npz(p, show=False, cmap="inferno")
        self.assertIsInstance(fig, plt.Figure)

    def test_string_path_accepted(self):
        p = make_snapshot(self.tmp / "snap.npz")
        fig, _ = plot_npz(str(p), show=False)
        self.assertIsInstance(fig, plt.Figure)


# =============================================================================
# make_gif_from_npz_files
# =============================================================================
class TestMakeGifFromNpzFiles(unittest.TestCase):

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)

    def tearDown(self):
        plt.close("all")
        self._tmp.cleanup()

    def test_creates_gif_file(self):
        paths = make_snapshot_dir(self.tmp / "run", n_steps=3)
        out = self.tmp / "out.gif"
        make_gif_from_npz_files(paths, out)
        self.assertTrue(out.exists())
        self.assertGreater(out.stat().st_size, 0)

    def test_returned_path_matches_gif_path(self):
        paths = make_snapshot_dir(self.tmp / "run", n_steps=2)
        out = self.tmp / "anim.gif"
        result = make_gif_from_npz_files(paths, out)
        self.assertEqual(result, out)

    def test_gif_is_valid_image(self):
        paths = make_snapshot_dir(self.tmp / "run", n_steps=3)
        out = self.tmp / "out.gif"
        make_gif_from_npz_files(paths, out)
        img = Image.open(out)
        self.assertEqual(img.format, "GIF")

    def test_gif_frame_count_matches_input(self):
        n = 4
        paths = make_snapshot_dir(self.tmp / "run", n_steps=n)
        out = self.tmp / "out.gif"
        make_gif_from_npz_files(paths, out)
        img = Image.open(out)
        frames = 0
        try:
            while True:
                frames += 1
                img.seek(img.tell() + 1)
        except EOFError:
            pass
        self.assertEqual(frames, n)

    def test_raises_for_empty_file_list(self):
        with self.assertRaises(ValueError):
            make_gif_from_npz_files([], self.tmp / "out.gif")

    def test_single_frame_gif_is_created(self):
        paths = [make_snapshot(self.tmp / "s0.npz")]
        out = self.tmp / "single.gif"
        make_gif_from_npz_files(paths, out)
        self.assertTrue(out.exists())

    def test_creates_parent_directory_if_needed(self):
        paths = make_snapshot_dir(self.tmp / "run", n_steps=2)
        out = self.tmp / "nested" / "deep" / "out.gif"
        make_gif_from_npz_files(paths, out)
        self.assertTrue(out.exists())

    def test_accepts_custom_fps(self):
        paths = make_snapshot_dir(self.tmp / "run", n_steps=2)
        out = self.tmp / "out.gif"
        make_gif_from_npz_files(paths, out, fps=12)
        self.assertTrue(out.exists())

    def test_accepts_custom_cmap(self):
        paths = make_snapshot_dir(self.tmp / "run", n_steps=2)
        out = self.tmp / "out.gif"
        make_gif_from_npz_files(paths, out, cmap="plasma")
        self.assertTrue(out.exists())

    def test_string_gif_path_accepted(self):
        paths = make_snapshot_dir(self.tmp / "run", n_steps=2)
        out = str(self.tmp / "out.gif")
        make_gif_from_npz_files(paths, out)
        self.assertTrue(Path(out).exists())


# =============================================================================
# _build_parser / CLI
# =============================================================================
class TestBuildParser(unittest.TestCase):

    def test_parser_is_returned(self):
        import argparse
        self.assertIsInstance(_build_parser(), argparse.ArgumentParser)

    def test_default_cmap_is_viridis(self):
        parser = _build_parser()
        args = parser.parse_args(["snap.npz"])
        self.assertEqual(args.cmap, "viridis")

    def test_default_fps_is_six(self):
        parser = _build_parser()
        args = parser.parse_args(["snap.npz"])
        self.assertEqual(args.fps, 6)

    def test_no_show_flag(self):
        parser = _build_parser()
        args = parser.parse_args(["snap.npz", "--no-show"])
        self.assertTrue(args.no_show)

    def test_save_flag(self):
        parser = _build_parser()
        args = parser.parse_args(["snap.npz", "--save", "out.png"])
        self.assertEqual(args.save_path, "out.png")

    def test_gif_flag(self):
        parser = _build_parser()
        args = parser.parse_args(["run/", "--gif", "out.gif"])
        self.assertEqual(args.gif_path, "out.gif")

    def test_cmap_flag(self):
        parser = _build_parser()
        args = parser.parse_args(["snap.npz", "--cmap", "plasma"])
        self.assertEqual(args.cmap, "plasma")

    def test_fps_flag(self):
        parser = _build_parser()
        args = parser.parse_args(["snap.npz", "--fps", "12"])
        self.assertEqual(args.fps, 12)


if __name__ == "__main__":
    unittest.main()

