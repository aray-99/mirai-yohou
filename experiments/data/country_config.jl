# 国別データ設定ローダ(Issue #3): experiments/data/countries/<ISO3>.toml を読む。
# 各データスクリプトから include される(M8 経由で二重 include されても安全なように
# 同値 const と関数再定義のみで構成)。
using TOML

const COUNTRY_TOML_DIR = joinpath(@__DIR__, "countries")

"countries/ に設定がある全 ISO3 コード(ソート順)"
list_countries() =
    sort([splitext(f)[1] for f in readdir(COUNTRY_TOML_DIR) if endswith(f, ".toml")])

"ISO3 の国別設定(Dict)。TOML の日付リテラルは Dates.Date で返る"
function load_country_config(iso3::AbstractString)
    path = joinpath(COUNTRY_TOML_DIR, "$(iso3).toml")
    isfile(path) || error("$path がありません。countries/README.md の書式で国別設定を作成してください")
    return TOML.parsefile(path)
end

"CLI 引数から ISO3 リストを取る(非フラグ引数。なければ countries/ の全国)"
country_args(args) = (c = [a for a in args if !startswith(a, "-")];
                      isempty(c) ? list_countries() : c)
