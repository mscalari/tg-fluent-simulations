"""Animação da variante "meio passo" da simulação MRF do conversor CT_V1.

A turbina (18 pás) é girada de 0° a 10° (meio passo, pois 360/18 = 20°) em
relação à bomba (20 pás) e ao estator (18 pás), em torno do eixo x. A vista é
axial (plano y-z, olhando ao longo de x), com as pás obtidas das faces de parede
da simulação SR000 (centro e vetor área de cada face). Não executa o Fluent.

Uso:
    PYTHONIOENCODING=utf-8 .venv241/Scripts/python.exe tools/animacao_meio_passo.py
    ... tools/animacao_meio_passo.py --limiar 0.85 --quadros 48 --fps 12
    ... tools/animacao_meio_passo.py --teste   (só o PNG de 0°)

Saídas (em results/conversor/CT_V1_MRF/graficos/):
    animacao_meio_passo.gif, animacao_meio_passo_0graus.png,
    animacao_meio_passo_10graus.png

Sentido positivo: regra da mão direita em torno de +x (y -> z):
    y' = y cos(t) - z sin(t);  z' = y sin(t) + z cos(t).
"""
import argparse
import io
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from PIL import Image

RAIZ = Path(__file__).resolve().parents[1]
DADOS = RAIZ / "results/conversor/CT_V1_MRF/SR000/dados"
SAIDA = RAIZ / "results/conversor/CT_V1_MRF/graficos"
BASE = "animacao_meio_passo"
TITULO = "Variante meio passo: turbina girada 10° em relação à bomba (MRF, rotor congelado)"
MAX_PONTOS = 60000
COR = {"bomba": "#1f77b4", "turbina": "#d62728", "estator": "#c8c8c8"}


def carregar_pa(peca, limiar):
    """Lê o CSV de parede de uma peça e retorna as faces de pá em (y, z).

    Args:
        peca (str): "bomba", "turbina" ou "estator".
        limiar (float): valor mínimo de |n·e_θ| para considerar a face de pá.

    Returns:
        numpy.ndarray: matriz (N, 2) com y e z (m) das faces classificadas como pá.

    Raises:
        FileNotFoundError: se o CSV da peça não existir.
    """
    arq = DADOS / f"SR000_parede_{peca}_area.csv"
    cols = ["y-coordinate", "z-coordinate", "face-area-magnitude",
            "y-face-area", "z-face-area"]
    d = pd.read_csv(arq, skipinitialspace=True, usecols=cols, dtype=np.float32)
    d.columns = d.columns.str.strip()
    y = d["y-coordinate"].to_numpy(np.float64)
    z = d["z-coordinate"].to_numpy(np.float64)
    r = np.hypot(y, z)
    mag = d["face-area-magnitude"].to_numpy(np.float64)
    ny = d["y-face-area"].to_numpy(np.float64) / mag
    nz = d["z-face-area"].to_numpy(np.float64) / mag
    comp_tang = np.abs(-ny * z + nz * y) / np.maximum(r, 1e-12)  # |n·e_θ|
    pa = comp_tang > limiar
    return np.column_stack([y[pa], z[pa]])


def subamostrar(pts, maximo=MAX_PONTOS):
    """Subamostra de forma determinística (passo fixo) até no máximo `maximo` pontos.

    Args:
        pts (numpy.ndarray): matriz (N, 2).
        maximo (int): número máximo de pontos.

    Returns:
        numpy.ndarray: matriz com no máximo `maximo` linhas.
    """
    if len(pts) <= maximo:
        return pts
    return pts[:: int(np.ceil(len(pts) / maximo))]


def rotacionar(pts, graus):
    """Gira pontos (y, z) em torno de +x pela regra da mão direita.

    Args:
        pts (numpy.ndarray): matriz (N, 2) com y e z.
        graus (float): ângulo de rotação em graus.

    Returns:
        numpy.ndarray: matriz (N, 2) com os pontos girados.
    """
    t = np.radians(graus)
    c, s = np.cos(t), np.sin(t)
    return np.column_stack([pts[:, 0] * c - pts[:, 1] * s,
                            pts[:, 0] * s + pts[:, 1] * c])


def centros_pas(pts, n_bins=720):
    """Estima nº de pás e ângulos de seus centros pelo espectro angular de uma coroa.

    O histograma angular das faces de pá é decomposto em harmônicos; o harmônico
    dominante (entre 10 e 30) dá o número de pás e sua fase dá a posição dos centros.

    Args:
        pts (numpy.ndarray): matriz (N, 2) de faces de pá.
        n_bins (int): número de divisões angulares do histograma.

    Returns:
        tuple: (n_pas, centros), com n_pas (int) e centros (numpy.ndarray) em graus.
    """
    r = np.hypot(pts[:, 0], pts[:, 1])
    sel = (r > 0.5 * r.max()) & (r < 0.9 * r.max())
    ang = np.degrees(np.arctan2(pts[sel, 1], pts[sel, 0])) % 360
    h, _ = np.histogram(ang, bins=n_bins, range=(0, 360))
    f = np.fft.rfft(h - h.mean())
    n = 10 + int(np.argmax(np.abs(f[10:31])))
    fase = -np.angle(f[n]) / n  # rad: posição do pico do harmônico n
    centros = (np.degrees(fase) + 360.0 * np.arange(n) / n) % 360
    return n, np.sort(centros)


def desenhar(pts, graus, lim, zoom):
    """Renderiza um quadro da animação e devolve a figura.

    Args:
        pts (dict): pontos por peça ("bomba", "turbina", "estator"), matrizes (N, 2).
        graus (float): rotação atual da turbina, em graus.
        lim (float): meia-largura (m) da vista completa.
        zoom (tuple): (centro_angulo_graus, r_min, r_max) do setor ampliado.

    Returns:
        matplotlib.figure.Figure: figura do quadro.
    """
    fig, axs = plt.subplots(1, 2, figsize=(11, 5.9), dpi=100)
    turb = rotacionar(pts["turbina"], graus)
    for ax in axs:
        ax.scatter(*pts["estator"].T, s=0.2, c=COR["estator"], lw=0)
        ax.scatter(*pts["bomba"].T, s=0.3, c=COR["bomba"], lw=0)
        ax.scatter(*turb.T, s=0.3, c=COR["turbina"], lw=0)
        ax.set_aspect("equal")
        ax.axis("off")
    axs[0].set_xlim(-lim, lim)
    axs[0].set_ylim(-lim, lim)
    axs[0].text(0, -lim * 1.02, f"rotação da turbina = {graus:4.1f}°",
                ha="center", va="top", fontsize=13, weight="bold")
    centro, rmin, rmax = zoom
    t = np.radians(centro)
    c = np.array([np.cos(t), np.sin(t)]) * (rmin + rmax) / 2
    meia = (rmax - rmin) * 0.62
    axs[1].set_xlim(c[0] - meia, c[0] + meia)
    axs[1].set_ylim(c[1] - meia, c[1] + meia)
    axs[1].set_title("zoom em setor de ~60°", fontsize=10)
    axs[0].set_title("vista axial (plano y-z, olhando ao longo de x)", fontsize=10)
    leg = [Line2D([], [], marker="o", ls="", c=COR["bomba"], label="bomba (fixa, 20 pás)"),
           Line2D([], [], marker="o", ls="", c=COR["turbina"],
                  label="turbina (gira 10° = meio passo)"),
           Line2D([], [], marker="o", ls="", c=COR["estator"], label="estator (fixo, 18 pás)")]
    fig.legend(handles=leg, loc="lower center", ncol=3, frameon=False, fontsize=10)
    fig.suptitle(TITULO, fontsize=12)
    fig.tight_layout(rect=(0, 0.06, 1, 0.95))
    return fig


def fig_para_imagem(fig):
    """Converte uma figura matplotlib em imagem PIL (RGB) e a fecha.

    Args:
        fig (matplotlib.figure.Figure): figura a converter.

    Returns:
        PIL.Image.Image: imagem RGB do quadro.
    """
    buf = io.BytesIO()
    fig.savefig(buf, format="png", facecolor="white")
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).convert("RGB")


def escolher_setor(centros_b, centros_t, rmax):
    """Escolhe o setor do zoom onde uma pá da turbina está mais alinhada a uma da bomba.

    Args:
        centros_b (numpy.ndarray): ângulos (graus) das pás da bomba.
        centros_t (numpy.ndarray): ângulos (graus) das pás da turbina.
        rmax (float): raio máximo (m) do conjunto.

    Returns:
        tuple: (ângulo central em graus, r_min, r_max) do setor.
    """
    dif = (centros_t[:, None] - centros_b[None, :] + 180) % 360 - 180
    i, j = np.unravel_index(np.argmin(np.abs(dif)), dif.shape)
    return (centros_t[i] + centros_b[j]) / 2, 0.45 * rmax, 0.95 * rmax


def principal(argv=None):
    """Gera o GIF e os PNGs do quadro inicial e final.

    Args:
        argv (list[str] | None): argumentos de linha de comando (None usa sys.argv).

    Returns:
        int: 0 em caso de sucesso.
    """
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--limiar", type=float, default=0.85, help="limiar de |n·e_θ|")
    ap.add_argument("--quadros", type=int, default=48, help="quadros do GIF (ida e volta)")
    ap.add_argument("--fps", type=float, default=12.0)
    ap.add_argument("--teste", action="store_true", help="só gera o PNG de 0°")
    a = ap.parse_args(argv)

    pts = {p: carregar_pa(p, a.limiar) for p in ("bomba", "turbina", "estator")}
    cent = {p: centros_pas(v) for p, v in pts.items()}
    for p, v in pts.items():
        print(f"{p}: {len(v)} faces de pá (limiar {a.limiar}), {cent[p][0]} pás detectadas")
    cb, ct = cent["bomba"][1], cent["turbina"][1]
    pts = {p: subamostrar(v) for p, v in pts.items()}
    rmax = max(np.hypot(*v.T).max() for v in pts.values())
    zoom = escolher_setor(cb, ct, rmax)
    lim = rmax * 1.05
    SAIDA.mkdir(parents=True, exist_ok=True)

    def quadro(g):
        """Renderiza o quadro para um ângulo.

        Args:
            g (float): ângulo da turbina em graus.

        Returns:
            PIL.Image.Image: quadro renderizado.
        """
        return fig_para_imagem(desenhar(pts, g, lim, zoom))

    quadro(0.0).save(SAIDA / f"{BASE}_0graus.png")
    if a.teste:
        return 0
    quadro(10.0).save(SAIDA / f"{BASE}_10graus.png")
    ida = np.linspace(0, 10, a.quadros // 2 + 1)
    angulos = np.r_[ida, ida[-2:0:-1]]
    imgs = [quadro(g).convert("P", palette=Image.ADAPTIVE, colors=64) for g in angulos]
    imgs[0].save(SAIDA / f"{BASE}.gif", save_all=True, append_images=imgs[1:],
                 duration=int(1000 / a.fps), loop=0, optimize=True)
    tam = (SAIDA / f"{BASE}.gif").stat().st_size / 1e6
    print(f"{len(imgs)} quadros; GIF {tam:.2f} MB")
    return 0


if __name__ == "__main__":
    raise SystemExit(principal())
