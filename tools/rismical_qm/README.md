# RISMiCal-QM

[![Language](https://img.shields.io/badge/Language-Python_3-blue.svg)](#)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](#)

[English](#english) | [日本語](#japanese)

<a id="english"></a>
## Overview
**RISMiCal-QM** is a highly optimized Python wrapper suite designed to perform robust QM/3D-RISM (KSDFT/3D-RISM, 3D-RISM-SCF, QM/MM/3D-RISM) calculations. It couples quantum chemistry engines (**Gaussian 16** or **PySCF**) with the 3D-RISM solver (**RISMiCal**).

## Key Features
* **Engine Flexibility:** Choose between `RISMiCal-QM-g16.py` (for Gaussian 16 users) and `RISMiCal-QM-PySCF.py` (for a modern, purely in-memory, I/O-free experience). Both wrappers share the exact same input format.
* **QM/MM/3D-RISM implemented:** By specifying `qmpart` in the input file, all atoms other than those specified will be treated as MM atoms. This enables QM/MM/3D-RISM calculations.
* **Rigorous Energy Decomposition:** Accurately separates the pure polarized QM internal energy ($E_{UU}$), the exact QM-MM interaction ($E_{QM-MM}$), and the QM/MM-Solvent interaction ($E_{QMMM-v}$) to calculate the exact Total Free Energy ($G_{TOTAL}$).
* **Excited State & Franck-Condon Support:** Seamlessly handles excited state potentials (TD-DFT). Supports Franck-Condon (FC) state calculations using frozen solvent environments via the `-FC` flag. All computed excitation energies are automatically extracted and summarized.
* **Automated Visualization Outputs:** Upon convergence, the script automatically generates `.cub` files for electrostatic potentials. In the PySCF wrapper, it also generates a single `.molden` file containing all molecular orbitals (MOs) for easy visualization in tools like Avogadro.


## Requirements
* **Python 3.6+**
* **NumPy**, **SciPy**
* **RISMiCal** (Command `rismical.x` must be in `$PATH`)
* For Gaussian version: **Gaussian 16** (`g16`, `formchk`, `cubegen` must be in `$PATH`)
* For PySCF version: **PySCF** (`pip install pyscf`)

## Usage
Run the script by passing the RISMiCal input file as an argument:


```
# For Gaussian 16 backend
python3 RISMiCal-QM-g16.py input.inp

# For PySCF backend
python3 RISMiCal-QM-PySCF.py input.inp

```

### Franck-Condon (FC) Mode

Use the `-FC` flag to perform a Franck-Condon (non-equilibrium) calculation:

```bash
python3 RISMiCal-QM-g16.py input.inp -FC

python3 RISMiCal-QM-PySCF.py input.inp -FC

```

---

## 概要

**RISMiCal-QM** は、QM/3D-RISM (KSDFT/3D-RISM, 3D-RISM-SCF, QM/MM/3D-RISM) 計算用Pythonラッパースイートです。
量子化学計算エンジン（**Gaussian 16** または **PySCF**）と3D-RISMソルバ（**RISMiCal**）を連成させ、溶質の電子状態、MM原子電荷、および溶媒の3D空間分布の間の自己無撞着場（SCF）サイクルを実行します。

## 主な特長

* **エンジンの選択:** 慣れ親しんだGaussian 16を利用する `RISMiCal-QM-g16.py` と、高速な完全インメモリ計算を実現する `RISMiCal-QM-PySCF.py` の2種類を提供。どちらも全く同じインプットフォーマットで動作します。
* **QM/MM/3D-RISMをサポート:** インプット中で`qmpart` を指定することで、指定された原子以外はMM原子として扱われます。これにより、QM/MM/3D-RISM計算が可能です。
* **厳密なエネルギー分解:** QM内部の分極エネルギー（$E_{UU}$）、QM-MM相互作用（$E_{QM-MM}$）、およびQM/MM-溶媒相互作用（$E_{QMMM-v}$）を厳密に分離して計算し、正しい全自由エネルギー（$G_{TOTAL}$）を出力します。
* **励起状態とFranck-Condon状態のサポート:** TD-DFT等を用いた平衡溶媒和計算にシームレスに対応。また `-FC` フラグを用いることで、凍結された溶媒場の中でのFranck-Condon状態（非平衡溶媒和）の計算を実行し、すべての励起エネルギーを自動集計します。
* **可視化ファイルの自動生成:** 収束後、静電ポテンシャルマップ（`.cub`）を自動生成します。PySCF版ではさらに、すべての分子軌道と軌道エネルギーをまとめた `.molden` ファイルを出力し、Avogadroなどでの直感的な可視化をサポートします。

## 動作環境

* **Python 3.6+**
* **NumPy**, **SciPy**
* **RISMiCal**（`rismical.x` にパスが通っていること）
* Gaussian版を利用する場合: **Gaussian 16**（`g16`, `formchk`, `cubegen` にパスが通っていること）
* PySCF版を利用する場合: **PySCF**（`pip install pyscf` でインストール）

## 使い方

RISMiCalのインプットファイルを引数に指定してスクリプトを実行します：

```bash
# Gaussian 16 を使用する場合
python3 RISMiCal-QM-g16.py input.inp

# PySCF を使用する場合
python3 RISMiCal-QM-PySCF.py input.inp

```

### Franck-Condon (FC) モード

`-FC` フラグを指定すると、Franck-Condon状態（非平衡溶媒和）の1点計算を実行します。

---

*Developed by Noriwo*


