# rismical-bind.py: 3D-RISM Binding Free Energy Analysis

## Overview
`rismical-bind.py` is a Python-based automation script for calculating ligand binding free energies using the **3D-RISM** (Three-Dimensional Reference Interaction Site Model) theory. It implements a workflow similar to **MM-PBSA**, combining molecular mechanics (MM) interaction energies with solvation free energies derived from RISMiCal.

## Features
- **Flexible Topology Support**: Handles AMBER (.prmtop) and GROMACS (.top) via ParmEd.
- **Advanced Trajectory Handling**: Supports various formats (.nc, .xtc, .mdcrd) via MDAnalysis.
- **Automated Workflow**: Iterates through snapshots, performs RMSD fitting/centering, generates RISMiCal inputs, executes calculations, and summarizes statistics.
- **Trajectory Export**: Can output aligned and centered trajectories for further analysis or visualization.
- **Robust Energy Calculation**: Computes intermolecular Coulomb and Lennard-Jones energies fast using Lorentz-Berthelot mixing rules via SciPy.

## Requirements
- Python 3.x
- **NumPy < 2.0** (Critical: NumPy 2.0+ is incompatible with current ParmEd versions)
- **ParmEd**
- **MDAnalysis**
- **SciPy**
- **RISMiCal** (`rismical.x` must be in your PATH)

## Usage
python rismical-bind.py <topology_file> <trajectory_file> <input_file>

Example:
python rismical-bind.py system.prmtop traj.nc rism.inp

## Input Configuration ($rismicalbind)
The script reads a custom namelist `$rismicalbind` within your RISMiCal `.inp` file:

| Parameter   | Description                                              | Default    |
| :---------- | :------------------------------------------------------- | :--------- |
| `host`      | Atom indices for the host (e.g., `1, -9200`).            | Required   |
| `guest`     | Atom indices for the guest (e.g., `9201, -9218`).        | 0 (None)   |
| `traj`      | Frame settings: `start, stop, interval`.                 | 1, last, 1 |
| `rmsfit`    | RMSD fitting target: `host`, `guest`, `complex`, `none`. | none       |
| `centering` | If `true`, centers the reference COM at (0,0,0).         | false      |
| `newtraj`   | If `true`, outputs a new trajectory file with fitted/centered coordinates (e.g., `traj.new.nc`). | false      |

## Formula
The binding free energy is calculated as the average over processed frames:
dG_bind = < E_MM + dSFE + dPC >

where:
- E_MM: Intermolecular MM interaction energy.
- dSFE: Solvation free energy change = SFE_SC(complex) - SFE_SC(host) - SFE_SC(guest)
- dPC: Change of Pressure correction term = Corr(complex) - Corr(host) - Corr(guest)

## Troubleshooting
- **UnicodeDecodeError**: If a binary NetCDF file has a `.mdcrd` extension, the script automatically attempts to force NetCDF reading to avoid MDAnalysis ASCII parsing errors.
- **ModuleNotFoundError (numpy.compat)**: Ensure you are using `numpy < 2.0`.

---

# rismical-bind.py: 3D-RISM 結合自由エネルギー解析

## 概要
`rismical-bind.py` は、**3D-RISM** (三次元参照相互作用部位モデル) 理論を用いてリガンドの結合自由エネルギーを計算するための Python 自動化スクリプトです。**MM-PBSA** 法と同様のワークフローを実装しており、分子力学 (MM) 相互作用エネルギーと、RISMiCal から得られる溶媒和自由エネルギーを組み合わせます。

## 主な機能
- **柔軟なトポロジ対応**: ParmEd を介して AMBER (.prmtop) や GROMACS (.top) を処理。
- **高度なトラジェクトリ操作**: MDAnalysis を用いて様々な形式 (.nc, .xtc, .mdcrd) をサポート。
- **自動化ワークフロー**: スナップショットの抽出、RMSD フィッティング、センタリング、RISMiCal 入力作成、実行、統計処理を一括で行います。
- **トラジェクトリ出力**: センタリング・フィッティング後の座標を新しいトラジェクトリとして保存可能です。
- **堅牢なエネルギー計算**: Lorentz-Berthelot 混合則に基づく、原子間のクーロンおよび Lennard-Jones 相互作用を高速に算出。

## 必要条件
- Python 3.x
- **NumPy < 2.0** (重要: NumPy 2.0 以降は現在の ParmEd と互換性がありません)
- **ParmEd**
- **MDAnalysis**
- **SciPy**
- **RISMiCal** (`rismical.x` が PATH に通っていること)

## 使用方法
python rismical-bind.py <トポロジファイル> <トラジェクトリファイル> <インプットファイル>

実行例:
python rismical-bind.py cb7-tbz_aq.top cb7-tbz_aq_production.nc cb7-tbz.inp

## インプット設定 ($rismicalbind)
RISMiCal の `.inp` ファイル内に独自のネームリスト `$rismicalbind` を記述します：

| パラメータ  | 説明                                                     | デフォルト |
| :---------- | :------------------------------------------------------- | :--------- |
| `host`      | ホスト分子の原子番号 (例: `1, -9200`)。                  | 指定必須   |
| `guest`     | ゲスト分子の原子番号 (例: `9201, -9218`)。               | 0 (なし)   |
| `traj`      | フレーム指定: `開始, 終了, 間隔`。                       | 1, last, 1 |
| `rmsfit`    | RMSD フィット対象: `host`, `guest`, `complex`, `none`。  | none       |
| `centering` | `true` の場合、リファレンスの重心を (0,0,0) に移動。     | false      |
| `newtraj`   | `true` の場合、補正後の座標を新しいトラジェクトリに出力（例: `sample.new.mdcrd`）。 | false      |

## 計算式
各フレームで算出された以下の値の平均を最終的な結合自由エネルギーとします：
dG_bind = < E_MM + dSFE + dPC >

内訳:
- E_MM: ホストーゲスト間の MM 相互作用エネルギー。
- dSFE: 溶媒和自由エネルギー変化 = SFE_SC(complex) - SFE_SC(host) - SFE_SC(guest)
- dPC: 圧力補正項の変化 = Corr(complex) - Corr(host) - Corr(guest)

## トラブルシューティング
- **UnicodeDecodeError**: バイナリの NetCDF トラジェクトリファイルが `.mdcrd` という拡張子を持っている場合、スクリプトは自動的に NetCDF リーダーを強制適用してエラーを回避します。
- **ModuleNotFoundError (numpy.compat)**: 環境の NumPy バージョンが 2.0 以上であることが原因です。`pip install "numpy<2.0"` を実行し、バージョンを下げてください。