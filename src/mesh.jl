struct fem_mesh{D<:AbstractDimension} <: AbstractMesh{D}
    nodes::Matrix{Float64}
    dofs::Int
    thick::Float64
end
function fem_mesh_gen(nodes; thick = nothing)
    dimension, node_number = size(nodes) # 2xnnodes或3xnnodes
    # 维度检查和厚度验证
    @assert !(dimension == 2 && isnothing(thick)) "For 2D problem, `thick` must be provided."
    thick = dimension == 2 ? thick : 1.0
    dofs = dimension * node_number
    D = dimension == 1 ? Dim1 :
        dimension == 2 ? Dim2 :
        dimension == 3 ? Dim3 :
        nothing
    fem_mesh{D}(nodes, dofs, thick)
end
