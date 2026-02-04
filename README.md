# Vector2Fold 🧬 — DNA Sequence to Protein 3D Structure Pipeline

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_svg.svg)](https://<あなたのアプリ名>.streamlit.app/)

## 🌟 プロジェクト概要
**Vector2Fold** は、ベクターなどの DNA 塩基配列を入力するだけで、タンパク質への翻訳、溶解性予測、そして最新の AI（ESMFold）を用いた 3D 構造予測までを一気通貫で行う Web アプリケーションです。

農学・バイオ研究において、実際にクローニングやタンパク質発現を行う前に、その配列が安定した立体構造を形成しうるかを数秒で視覚化し、実験の試行錯誤を減らすことを目的としています。

[Image showing a DNA sequence being translated into a 3D protein structure]

## 🧪 主な機能
- **DNA to Protein 翻訳**: 遺伝子配列（DNA）をアミノ酸配列へ自動変換。
- **溶解性・安定性指標の提示**: ESMFold による信頼度スコア（pLDDT）に基づき、構造の安定性を評価。
- **3D インタラクティブ表示**: `py3Dmol` を用いた立体構造の可視化。pLDDT スコアに応じた色分け（青：高信頼度、赤：低信頼度）に対応。
- **正規化ロジック搭載**: API レスポンスのスケールを自動判別し、常に 100% 満点で評価を表示。

## 🛠 テクニカルスタック
- **Frontend**: Streamlit
- **Bioinformatics**: Biopython
- **AI Model**: ESMFold (ESM-2) API
- **Visualization**: py3Dmol / stmol
- **Language**: Python 3.12

## 🚀 使い方
1. [デプロイ済みアプリ](https://<あなたのアプリ名>.streamlit.app/)にアクセスします。
2. DNA 配列（例：ユビキチンやインスリンの CDS 部分）を入力欄にペーストします。
3. 「解析開始」ボタンをクリックします。
4. 数秒待つと、翻訳されたアミノ酸配列、予測スコア、および 3D 構造が表示されます。

## 📂 リポジトリ構成
- `main.py`: アプリケーションのメインスクリプト
- `requirements.txt`: 実行に必要なライブラリ一覧（Cython 含む）
- `LICENSE`: MIT License

## 👤 作成者
- **氏名**: [あなたの名前]
- **専門**: 博士（農学）
- **研究背景**: 魚類免疫学、次世代シーケンサー解析、抗体医薬品 R&D
- **関心分野**: バイオインフォマティクス / Python による解析ツール開発

## 📄 ライセンス
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.