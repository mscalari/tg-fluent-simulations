"""Mapas de campo do conversor CT_V1 (MRF) para a Seção 4.5 do TG, a partir das imagens já exportadas pelo Fluent.

As imagens do Fluent (``results/conversor/CT_V1_MRF/SRxxx/graficos/*.png``, copiadas para ``C_TG/img`` do
Obsidian com o prefixo ``CT_V1_MRF_``) têm a legenda em inglês, o nome da superfície ("c_mer", "pl_bt") e muito
espaço em branco. Este script:

* recorta a região do campo (toda a área não branca à direita da legenda; nos planos meridionais, só a metade
  superior, y > 0);
* recorta a faixa de cores da legenda original e a redesenha como barra de cores do matplotlib, com o título em
  português, a unidade e os valores mínimo e máximo lidos da legenda do Fluent (escala linear, como no Fluent);
* monta os painéis lado a lado e grava em ``tg-fem-unicamp/Figuras/conversor_campo_*.png``.

Os extremos de cada legenda (``ESCALAS``) foram lidos das próprias imagens; conferir se uma imagem for refeita.
O campo não é reamostrado: os pixels são os do Fluent.

Uso (rodar com ``PYTHONIOENCODING=utf-8`` e o ``.venv241``)::

    python tools/figuras_tg_campos_conversor.py              gera todas as figuras
    python tools/figuras_tg_campos_conversor.py vtheta pas   gera só as escolhidas (nomes em ``FIGURAS``)
"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MaxNLocator
from PIL import Image

RAIZ = Path(__file__).resolve().parents[1]
IMG = Path("C:/Users/murilo/Repositorios/obsidian-vault-murilo-scalari/CLAUDE/C_TG/img")
SAIDA = RAIZ.parent / "tg-fem-unicamp" / "Figuras"
SENS = RAIZ / "results" / "conversor" / "CT_V1_MRF" / "sens_campos"   # mapas de campo_sensibilidade_conversor.py

# (vmin, vmax) lidos da legenda de cada imagem do Fluent, nas unidades do Fluent (Pa, m/s, adimensional)
ESCALAS = {
    "SR000_01_meridional_pressao": (-3.26e6, 1.38e6), "SR090_01_meridional_pressao": (-2.38e6, 6.82e5),
    "SR000_05_meridional_v_axial": (-92.3, 37.4), "SR090_05_meridional_v_axial": (-57.4, 16.2),
    "SR000_07_meridional_v_tangencial": (-41.0, 91.5), "SR040_07_meridional_v_tangencial": (-15.6, 99.3),
    "SR090_07_meridional_v_tangencial": (-13.3, 79.2), "SR100_07_meridional_v_tangencial": (-16.8, 78.3),
    "SR000_09_meridional_mut_mu": (0.836, 2.50e3), "SR090_09_meridional_mut_mu": (1.0, 445.0),
    "SR000_4_bomba_pressao_de_mais_x": (-5.41e6, 3.65e6), "SR090_4_bomba_pressao_de_mais_x": (-3.14e6, 7.63e5),
    "SR000_4_turbina_pressao_de_menos_x": (-7.72e6, 4.00e6), "SR090_4_turbina_pressao_de_menos_x": (-1.62e6, 7.18e5),
    "SR000_6_pathlines_meridional": (0.0, 134.0), "SR090_6_pathlines_meridional": (0.0, 86.2),
    "SR000_5_if_bomba_turbina_lado_bomba_v_axial": (-56.5, 58.2),
    "SR090_5_if_bomba_turbina_lado_bomba_v_axial": (-22.2, 21.0),
    "SR000_4_bomba_yplus_de_mais_x": (2.11, 1160.0), "SR000_4_turbina_yplus_de_menos_x": (1.28, 1140.0),
}

plt.rcParams.update({"font.size": 8, "savefig.dpi": 250, "mathtext.fontset": "cm"})


def fmt_virgula(v, _):
    """Formata um valor de marcação com vírgula decimal e sinal de menos tipográfico.

    Args:
        v (float): Valor.
        _ (int): Posição (não usada).

    Returns:
        str: Texto do rótulo.
    """
    return f"{v:.4g}".replace(".", ",").replace("-", "−")


def ler(nome):
    """Lê uma imagem do Fluent como array RGB.

    Args:
        nome (str): Nome sem o prefixo ``CT_V1_MRF_`` e sem extensão; ``"sens:<arquivo>"`` lê de ``SENS``.

    Returns:
        np.ndarray: Imagem (altura, largura, 3), uint8.

    Raises:
        FileNotFoundError: Se a imagem não existir.
    """
    arq = SENS / f"{nome[5:]}.png" if nome.startswith("sens:") else IMG / f"CT_V1_MRF_{nome}.png"
    return np.array(Image.open(arq).convert("RGB"))


def escala_de(nome):
    """Extremos da legenda de uma imagem: ``ESCALAS`` ou, nas de sensibilidade, a faixa fixa de ``faixas.json``.

    Args:
        nome (str): Nome da imagem (``"sens:<caso>_<campo>"`` para os mapas de sensibilidade).

    Returns:
        tuple[float, float]: (mínimo, máximo) nas unidades do Fluent.
    """
    if nome.startswith("sens:"):
        faixas = json.loads((SENS / "faixas.json").read_text())["faixas"]
        campo = next(c for c in sorted(faixas, key=len, reverse=True) if nome.endswith("_" + c))
        return tuple(faixas[campo])
    return ESCALAS[nome]


def faixa_cores(a):
    """Localiza e recorta a faixa de cores da legenda (pixels saturados nas 300 colunas da esquerda).

    Args:
        a (np.ndarray): Imagem RGB.

    Returns:
        np.ndarray: Faixa (altura, largura, 3), com o máximo da escala no topo.
    """
    e = a[:, :300].astype(int)
    ys, xs = np.where((e.max(2) - e.min(2)) > 80)
    return a[ys.min() + 1:ys.max(), xs.min() + 4:xs.max() - 4]


def recorte_campo(a, modo="todo", margem=6):
    """Recorta a região do campo, à direita da legenda.

    Args:
        a (np.ndarray): Imagem RGB.
        modo (str): ``"todo"`` (caixa de todos os pixels não brancos) ou ``"superior"`` (só o bloco de cima,
            separado do de baixo por linhas brancas; nos planos meridionais é a metade y > 0).
        margem (int): Pixels de folga em volta.

    Returns:
        np.ndarray: Recorte RGB.
    """
    b = a[:, 320:]
    nb = (b < 245).any(2)
    linhas = np.where(nb.any(1))[0]
    y0, y1 = linhas.min(), linhas.max()
    if modo == "superior":
        vazias = [y for y in range(y0, y1) if not nb[y].any()]
        if vazias:
            y1 = vazias[0] - 1
    cols = np.where(nb[y0:y1 + 1].any(0))[0]
    x0, x1 = cols.min(), cols.max()
    return b[max(y0 - margem, 0):y1 + margem, max(x0 - margem, 0):x1 + margem]


def painel(fig, spec, nome, titulo, rotulo, escala=1.0, modo="todo", textos=()):
    """Desenha um mapa recortado com a barra de cores refeita ao lado.

    Args:
        fig (matplotlib.figure.Figure): Figura.
        spec (matplotlib.gridspec.SubplotSpec): Espaço do painel.
        nome (str): Imagem do Fluent (chave de ``ESCALAS``).
        titulo (str): Título do painel.
        rotulo (str): Rótulo da barra de cores (grandeza e unidade).
        escala (float): Fator aplicado aos extremos (por exemplo 1e-6 para Pa -> MPa).
        modo (str): Modo de recorte (``recorte_campo``).
        textos (tuple): Anotações ``(x, y, texto, ha)`` em fração do eixo do mapa.
    """
    a = ler(nome)
    vmin, vmax = (v * escala for v in escala_de(nome))
    g = spec.subgridspec(1, 2, width_ratios=[1, 0.055], wspace=0.04)
    ax = fig.add_subplot(g[0])
    ax.imshow(recorte_campo(a, modo))
    ax.set_axis_off()
    ax.set_title(titulo, fontsize=8.5)
    for x, y, t, ha in textos:
        ax.text(x, y, t, transform=ax.transAxes, ha=ha, va="top", fontsize=7, color="0.15")
    cb = fig.add_subplot(g[1])
    cb.imshow(faixa_cores(a), extent=[0, 1, vmin, vmax], aspect="auto")
    cb.set_xticks([])
    cb.yaxis.tick_right()
    cb.yaxis.set_label_position("right")
    cb.yaxis.set_major_locator(MaxNLocator(nbins=5))
    cb.set_ylim(vmin, vmax)
    cb.yaxis.set_major_formatter(FuncFormatter(fmt_virgula))
    cb.tick_params(labelsize=6.5, length=2)
    cb.set_ylabel(rotulo, fontsize=7.5)


MERID = ((0.02, 0.99, "bomba", "left"), (0.98, 0.99, "turbina", "right"))


def salvar(fig, nome):
    """Grava a figura em ``SAIDA``.

    Args:
        fig (matplotlib.figure.Figure): Figura.
        nome (str): Nome sem extensão.
    """
    SAIDA.mkdir(parents=True, exist_ok=True)
    fig.savefig(SAIDA / f"{nome}.png", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("gravado", SAIDA / f"{nome}.png")


def fig_meridional():
    """Pressão estática e velocidade axial no plano meridional, em SR 0 e 0,9."""
    fig = plt.figure(figsize=(6.6, 5.6))
    gs = fig.add_gridspec(2, 2, hspace=0.18, wspace=0.28)
    painel(fig, gs[0, 0], "SR000_01_meridional_pressao", "(a) Pressão estática, SR = 0", "$p$ [MPa]", 1e-6,
           "superior", MERID)
    painel(fig, gs[0, 1], "SR090_01_meridional_pressao", "(b) Pressão estática, SR = 0,9", "$p$ [MPa]", 1e-6,
           "superior", MERID)
    painel(fig, gs[1, 0], "SR000_05_meridional_v_axial", "(c) Velocidade axial, SR = 0", "$c_x$ [m/s]", 1.0,
           "superior", MERID)
    painel(fig, gs[1, 1], "SR090_05_meridional_v_axial", "(d) Velocidade axial, SR = 0,9", "$c_x$ [m/s]", 1.0,
           "superior", MERID)
    salvar(fig, "conversor_campo_meridional")


def fig_vtheta():
    """Velocidade tangencial absoluta no plano meridional para SR 0, 0,4, 0,9 e 1,0."""
    fig = plt.figure(figsize=(6.6, 5.6))
    gs = fig.add_gridspec(2, 2, hspace=0.18, wspace=0.28)
    for k, (sr, rot) in enumerate((("SR000", "0"), ("SR040", "0,4"), ("SR090", "0,9"), ("SR100", "1,0"))):
        painel(fig, gs[k // 2, k % 2], f"{sr}_07_meridional_v_tangencial", f"({'abcd'[k]}) SR = {rot}",
               r"$c_\theta$ [m/s]", 1.0, "superior", MERID)
    salvar(fig, "conversor_campo_vtheta")


def fig_mut():
    """Razão de viscosidade turbulenta no plano meridional, em SR 0 e 0,9."""
    fig = plt.figure(figsize=(6.6, 2.9))
    gs = fig.add_gridspec(1, 2, wspace=0.28)
    painel(fig, gs[0], "SR000_09_meridional_mut_mu", "(a) SR = 0", r"$\mu_t/\mu$", 1.0, "superior", MERID)
    painel(fig, gs[1], "SR090_09_meridional_mut_mu", "(b) SR = 0,9", r"$\mu_t/\mu$", 1.0, "superior", MERID)
    salvar(fig, "conversor_campo_mut")


def fig_pas():
    """Pressão estática nas paredes da bomba (vista do lado da turbina) e da turbina (vista do lado da bomba)."""
    fig = plt.figure(figsize=(6.6, 5.8))
    gs = fig.add_gridspec(2, 2, hspace=0.15, wspace=0.28)
    painel(fig, gs[0, 0], "SR000_4_bomba_pressao_de_mais_x", "(a) Bomba, SR = 0", "$p$ [MPa]", 1e-6)
    painel(fig, gs[0, 1], "SR090_4_bomba_pressao_de_mais_x", "(b) Bomba, SR = 0,9", "$p$ [MPa]", 1e-6)
    painel(fig, gs[1, 0], "SR000_4_turbina_pressao_de_menos_x", "(c) Turbina, SR = 0", "$p$ [MPa]", 1e-6)
    painel(fig, gs[1, 1], "SR090_4_turbina_pressao_de_menos_x", "(d) Turbina, SR = 0,9", "$p$ [MPa]", 1e-6)
    salvar(fig, "conversor_campo_pas")


def fig_trajetorias():
    """Linhas de trajetória lançadas da interface bomba-turbina, vistas perpendicularmente ao eixo, SR 0 e 0,9."""
    fig = plt.figure(figsize=(6.0, 4.4))
    gs = fig.add_gridspec(1, 2, wspace=0.3)
    txt = ((0.02, 0.99, "bomba", "left"), (0.98, 0.99, "turbina", "right"))
    painel(fig, gs[0], "SR000_6_pathlines_meridional", "(a) SR = 0", "$|c|$ [m/s]", 1.0, "todo", txt)
    painel(fig, gs[1], "SR090_6_pathlines_meridional", "(b) SR = 0,9", "$|c|$ [m/s]", 1.0, "todo", txt)
    salvar(fig, "conversor_campo_trajetorias")


def fig_if_bt():
    """Velocidade axial no lado da bomba da interface bomba-turbina, SR 0 e 0,9 (positivo: da bomba para a turbina)."""
    fig = plt.figure(figsize=(6.6, 3.0))
    gs = fig.add_gridspec(1, 2, wspace=0.3)
    painel(fig, gs[0], "SR000_5_if_bomba_turbina_lado_bomba_v_axial", "(a) SR = 0", "$c_x$ [m/s]")
    painel(fig, gs[1], "SR090_5_if_bomba_turbina_lado_bomba_v_axial", "(b) SR = 0,9", "$c_x$ [m/s]")
    salvar(fig, "conversor_campo_if_bt")


def fig_yplus():
    """y⁺ nas paredes da bomba e da turbina no estol."""
    fig = plt.figure(figsize=(6.6, 3.0))
    gs = fig.add_gridspec(1, 2, wspace=0.3)
    painel(fig, gs[0], "SR000_4_bomba_yplus_de_mais_x", "(a) Bomba", "$y^+$")
    painel(fig, gs[1], "SR000_4_turbina_yplus_de_menos_x", "(b) Turbina", "$y^+$")
    salvar(fig, "conversor_campo_yplus")


def fig_sens_mut():
    """Razão de viscosidade turbulenta no plano meridional, no estol: modelos de turbulência e esquema de 1ª ordem."""
    fig = plt.figure(figsize=(6.6, 5.6))
    gs = fig.add_gridspec(2, 2, hspace=0.18, wspace=0.28)
    for k, (caso, tit) in enumerate((("base", "(a) SST, segunda ordem (base)"),
                                     ("sst_sem_cc", "(b) SST sem correção de curvatura"),
                                     ("ke_realizable", r"(c) $k$-$\varepsilon$ Realizable"),
                                     ("fou", "(d) SST, primeira ordem"))):
        painel(fig, gs[k // 2, k % 2], f"sens:{caso}_mut", tit, r"$\mu_t/\mu$", 1.0, "superior", MERID)
    salvar(fig, "conversor_campo_sens_mut")


def fig_sens_vtheta():
    """Velocidade tangencial no plano meridional, no estol: SOU × FOU na G0 e SIMPLE × Coupled na G2."""
    fig = plt.figure(figsize=(6.6, 5.6))
    gs = fig.add_gridspec(2, 2, hspace=0.18, wspace=0.28)
    for k, (caso, tit) in enumerate((("base", "(a) G0, segunda ordem (base)"), ("fou", "(b) G0, primeira ordem"),
                                     ("g2_simple", "(c) G2, SIMPLE"), ("g2_coupled", "(d) G2, Coupled"))):
        painel(fig, gs[k // 2, k % 2], f"sens:{caso}_vtheta", tit, r"$c_\theta$ [m/s]", 1.0, "superior", MERID)
    salvar(fig, "conversor_campo_sens_vtheta")


def fig_sens_ifbt():
    """Velocidade axial no lado da bomba da interface bomba-turbina, no estol, na posição original e no meio passo."""
    fig = plt.figure(figsize=(6.6, 3.0))
    gs = fig.add_gridspec(1, 2, wspace=0.3)
    painel(fig, gs[0], "sens:base_ifbt_vx", "(a) Posição original (base)", "$c_x$ [m/s]")
    painel(fig, gs[1], "sens:meio_passo_ifbt_vx", "(b) Turbina girada de meio passo", "$c_x$ [m/s]")
    salvar(fig, "conversor_campo_sens_ifbt")


FIGURAS = {"merid": fig_meridional, "vtheta": fig_vtheta, "mut": fig_mut, "pas": fig_pas,
           "traj": fig_trajetorias, "ifbt": fig_if_bt, "yplus": fig_yplus,
           "sens_mut": fig_sens_mut, "sens_vtheta": fig_sens_vtheta, "sens_ifbt": fig_sens_ifbt}

if __name__ == "__main__":
    for nome in (sys.argv[1:] or FIGURAS):
        FIGURAS[nome]()
