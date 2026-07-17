---
name: 作業タスク(自己完結 Issue)
about: 未来予報の委譲可能な作業単位。この Issue の記述だけで作業が完了できること
labels: sonnet-ok
---

## 目的

<!-- なぜやるか。1〜3 文。関連する DECISIONS 番号・設計文書の節を挙げる -->

## 対象

<!-- 触るファイル・作らないファイル。読むべき仕様の箇所(docs/SPEC.md §X など) -->

## 完了条件

<!-- チェックリスト形式。テスト green(julia --project -e 'using Pkg; Pkg.test()')は常に含める -->

- [ ]

## 報告形式

<!-- PR に何を書くか。数値結果が出る場合はその形式 -->

## 環境要件

<!-- ACLED 認証の要否(キャッシュ済み CSV で足りるか)、長時間ランの有無と目安時間 -->

## 制約(全 Issue 共通)

- 凍結済み合格基準・設定値(acceptance_thresholds.toml、M8_frozen_config.toml、DECISIONS で凍結された値)は変更しない
- Conventional Commits(件名英語、scope は CLAUDE.md 参照)。DECISIONS 参照は本文末尾の `Refs:` 行
- テスト red でコミットしない。実験成果物には来歴サイドカーを付す
