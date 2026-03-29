import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

DEFAULT_FIELDS = ("phi", "psi", "ngf", "mf", "ml", "mtb")


def load_simulation_npz(npz_path):
    """Load a simulation .npz file produced by main.py.

    Returns a dict keyed by variable name.
    """
    npz_path = Path(npz_path)
    if not npz_path.exists():
        raise FileNotFoundError(f"File not found: {npz_path}")

    with np.load(npz_path, allow_pickle=False) as data:
        return {key: data[key] for key in data.files}


def _step_sort_key(path_obj):
    """Sort paths by the trailing number in the filename when available."""
    matches = re.findall(r"(\d+)", path_obj.stem)
    if matches:
        return int(matches[-1]), path_obj.name
    return float("inf"), path_obj.name


def collect_npz_files(source):
    """Collect .npz files from a file, directory, or glob pattern."""
    source_str = str(source)
    path = Path(source_str)

    if path.is_file():
        return [path]

    if path.is_dir():
        files = sorted(path.glob("*.npz"), key=_step_sort_key)
        if not files:
            raise ValueError(f"No .npz files found in directory: {path}")
        return files

    # Fallback: treat source as glob pattern
    parent = path.parent if str(path.parent) not in ("", ".") else Path(".")
    files = sorted(parent.glob(path.name), key=_step_sort_key)
    if not files:
        raise ValueError(f"No files matched pattern: {source}")
    return files


def plot_simulation_snapshot(data, title=None, cmap="viridis"):
    """Plot the fields saved by main.py in a compact figure.

    Parameters
    ----------
    data : dict
        Dictionary returned by ``load_simulation_npz``.
    title : str | None
        Optional figure title.
    cmap : str
        Matplotlib colormap used for scalar fields.

    Returns
    -------
    (fig, axes)
        Matplotlib figure and axes.
    """
    fig, axes = plt.subplots(2, 4, figsize=(16, 8), constrained_layout=True)
    axes = axes.ravel()

    used_axes = 0
    for field in DEFAULT_FIELDS:
        if field in data:
            ax = axes[used_axes]
            im = ax.imshow(data[field], origin="lower", cmap=cmap)
            ax.set_title(field)
            ax.set_xticks([])
            ax.set_yticks([])
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            used_axes += 1

    if "v_m" in data:
        ax = axes[used_axes]
        vm = data["v_m"]
        if vm.ndim == 3 and vm.shape[0] == 2:
            speed = np.sqrt(vm[0] ** 2 + vm[1] ** 2)
            im = ax.imshow(speed, origin="lower", cmap="magma")
            ax.set_title("v_m magnitude")
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        else:
            im = ax.imshow(np.asarray(vm), origin="lower", cmap="magma")
            ax.set_title("v_m")
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        ax.set_xticks([])
        ax.set_yticks([])
        used_axes += 1

    for i in range(used_axes, len(axes)):
        axes[i].axis("off")

    if title:
        fig.suptitle(title)

    return fig, axes


def plot_npz(npz_path, save_path=None, show=True, cmap="viridis"):
    """Convenience wrapper: load a .npz and plot all relevant fields."""
    data = load_simulation_npz(npz_path)
    fig, axes = plot_simulation_snapshot(data, title=Path(npz_path).name, cmap=cmap)

    if save_path is not None:
        fig.savefig(save_path, dpi=200)

    if show:
        plt.show()

    return fig, axes


def make_gif_from_npz_files(npz_files, gif_path, fps=6, cmap="viridis"):
    """Create a GIF from an ordered list of NPZ simulation snapshots."""
    if not npz_files:
        raise ValueError("No NPZ files were provided.")

    frames = []
    frame_duration_ms = max(1, int(1000 / max(1, fps)))

    for npz_file in npz_files:
        data = load_simulation_npz(npz_file)
        fig, _ = plot_simulation_snapshot(data, title=Path(npz_file).name, cmap=cmap)

        fig.canvas.draw()
        rgba = np.asarray(fig.canvas.buffer_rgba())

        # Backends may return RGBA as (H, W, 4), packed (H, W), or flat bytes.
        if rgba.ndim == 3 and rgba.shape[-1] == 4:
            rgb = rgba[:, :, :3]
        elif rgba.ndim == 2:
            packed = np.asarray(rgba, dtype=np.uint32)
            rgb = packed.view(np.uint8).reshape(packed.shape[0], packed.shape[1], 4)[:, :, :3]
        elif rgba.ndim == 1:
            width, height = fig.canvas.get_width_height()
            total_pixels = rgba.size // 4
            base_pixels = width * height
            scale = int(round((total_pixels / base_pixels) ** 0.5)) if base_pixels else 1
            scale = max(1, scale)
            rgb = rgba.reshape(height * scale, width * scale, 4)[:, :, :3]
        else:
            raise ValueError(f"Unsupported canvas buffer shape: {rgba.shape}")

        frames.append(Image.fromarray(np.ascontiguousarray(rgb)))
        plt.close(fig)

    gif_path = Path(gif_path)
    gif_path.parent.mkdir(parents=True, exist_ok=True)
    frames[0].save(
        gif_path,
        save_all=True,
        append_images=frames[1:],
        duration=frame_duration_ms,
        loop=0,
    )
    return gif_path


def _build_parser():
    parser = argparse.ArgumentParser(
        description="Plot simulation .npz files generated by pyaxon/main.py"
    )
    parser.add_argument(
        "npz_path",
        nargs="?",
        help="Single .npz file for plotting, or directory/glob source when using --gif",
    )
    parser.add_argument("--save", dest="save_path", default=None, help="Optional output image path")
    parser.add_argument("--no-show", action="store_true", help="Do not open an interactive window")
    parser.add_argument("--cmap", default="viridis", help="Colormap for scalar fields")
    parser.add_argument("--gif", dest="gif_path", default=None, help="Output GIF path")
    parser.add_argument("--fps", type=int, default=6, help="GIF frame rate")
    return parser


def main():
    parser = _build_parser()
    args = parser.parse_args()

    if args.gif_path:
        if args.npz_path is None:
            parser.error("npz_path is required when using --gif (file, directory, or glob pattern)")
        npz_files = collect_npz_files(args.npz_path)
        out = make_gif_from_npz_files(npz_files, args.gif_path, fps=args.fps, cmap=args.cmap)
        print(f"Saved GIF: {out}")
        return

    if args.npz_path is None:
        parser.error("npz_path is required for plotting")

    plot_npz(args.npz_path, save_path=args.save_path, show=not args.no_show, cmap=args.cmap)


if __name__ == "__main__":
    main()
