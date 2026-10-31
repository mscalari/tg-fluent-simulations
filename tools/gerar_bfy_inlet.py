"""Gera ``refs/bfy_inlet.prof``: lei de potência 1/7 nas duas paredes do canal de entrada do degrau turbulento.

Uso: python tools/gerar_bfy_inlet.py
"""
import numpy as np
from casos import ROOT

DELTA, Y0, Y1, N = 1.4, 1.0, 9.0, 801


def perfil(y, delta=DELTA):
    """Velocidade da lei 1/7 com núcleo U = 1, tomando a distância à parede mais próxima.

    Args:
        y (np.ndarray): Coordenada [m], em [1, 9].
        delta (float, optional): Espessura da camada-limite [m]. Padrão 1,4.

    Returns:
        np.ndarray: u(y) [m/s].
    """
    d = np.maximum(np.minimum(y - Y0, Y1 - y), 0.0)
    return np.minimum(d / delta, 1.0) ** (1 / 7)


def bloco(nome, v):
    """Formata uma coluna do perfil Fluent.

    Args:
        nome (str): Nome do campo.
        v (np.ndarray): Valores.

    Returns:
        str: Texto da coluna.
    """
    lin = [" ".join(f"{a:.9E}" for a in v[i:i + 5]) for i in range(0, len(v), 5)]
    return f" ({nome}\n   " + "\n   ".join(lin) + "\n )\n"


def main():
    """Grava o perfil com pontos concentrados nas paredes (malha em cosseno) e confere os limites."""
    y = Y0 + (Y1 - Y0) * (1 - np.cos(np.linspace(0, np.pi, N))) / 2
    u = perfil(y)
    assert u[0] == 0 and abs(u[-1]) < 1e-1 and abs(u[N // 2] - 1) < 1e-12
    assert abs(perfil(np.array([2.4]))[0] - 1) < 1e-12
    txt = f"((bfy-inlet point {N})\n" + bloco("x", np.full(N, -10.0)) + bloco("y", y) + bloco("u", u) + ")\n"
    (ROOT / "refs" / "bfy_inlet.prof").write_text(txt)
    print("ū =", np.trapezoid(u, y) / (Y1 - Y0))


if __name__ == "__main__":
    main()
