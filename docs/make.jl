using Documenter, FEM
DocMeta.setdocmeta!(FEM, :DocTestSetup, :(using FEM); recursive=true)
makedocs(
    modules = [FEM],
    sitename = "FEM.jl",
    authors  = "Sun Haohang",
    pages = [
        "index.md",
        "代码文档" => [
            "形函数" => joinpath("man", "shape.md"),
            "其他" => joinpath("man", "other.md")
        ],
        "理论推导备忘" => [
            "弹性理论及有限元" => joinpath("theory", "elasticity.md"),
            "梯度损伤模型" => joinpath("theory", "gradient.md")
            ]
    ]
)
