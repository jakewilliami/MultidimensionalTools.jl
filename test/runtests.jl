using MultidimensionalTools
using Test

function aresamearrays(A::AbstractVector{T}, B::AbstractVector{R}) where {T, R}
	A === B && return true
	length(A) ≠ length(B) && return false
	
	for (a, b) in zip(sort(A), sort(B))
		if a ≠ b
			return false
		end
	end
	
	return true
end
aresamearrays(A::AbstractArray{T, N}, B::AbstractArray{R, M}) where {T, R, N, M} =
	aresamearrays(reshape(A, :), reshape(B, :))

# TESTS

@time @testset "MultidimensionalTools.jl Types" begin
	# TYPES
	@test (0, 0, 0) isa AbstractIndex{3}
	@test CartesianIndex(1, 2, 3, 4, 5) isa AbstractIndex{5}
	@test [(0, 0), (1, 1), (2, 2)] isa AbstractIndices{2}
	@test CartesianIndex.([(3, 6, 9, 12), (4, 4, 1, 0)]) isa AbstractIndices{4}
	@test CartesianIndices([1, 2, 3, 4]) isa AbstractIndices{1}
	@test CartesianIndices(rand(Int8, 3, 3, 3, 3, 3, 3, 3)) isa AbstractIndices{7}
	@test (0, 0, 0) isa AbstractIndexOrIndices{3}
	@test CartesianIndex(1, 2, 3, 4, 5) isa AbstractIndexOrIndices{5}
	@test [(0, 0), (1, 1), (2, 2)] isa AbstractIndexOrIndices{2}
	@test CartesianIndex.([(3, 6, 9), (4, 4, 1)]) isa AbstractIndexOrIndices{3}
	@test CartesianIndices([1, 2, 3, 4]) isa AbstractIndexOrIndices{1}
	@test CartesianIndices(rand(Int8, 3, 3, 3, 3, 3, 3, 3)) isa AbstractIndexOrIndices{7}
end

@time @testset "MultidimensionalTools.jl Main" begin
	# MAIN
    @test n_adjacencies(2) == 8
    @test n_adjacencies(3) == 26
    @test n_adjacencies(17) == 129140162
	@test n_adjacencies(cat([1 1; 1 1], [1 1; 1 1], [1 1; 1 1], dims = 3)) == n_adjacencies(3)
    @test aresamearrays(get_directions(3), [(-1, -1, -1), (0, -1, -1), (1, -1, -1), (-1, 0, -1), (0, 0, -1), (1, 0, -1), (-1, 1, -1), (0, 1, -1), (1, 1, -1), (-1, -1, 0), (0, -1, 0), (1, -1, 0), (-1, 0, 0), (1, 0, 0), (-1, 1, 0), (0, 1, 0), (1, 1, 0), (-1, -1, 1), (0, -1, 1), (1, -1, 1), (-1, 0, 1), (0, 0, 1), (1, 0, 1), (-1, 1, 1), (0, 1, 1), (1, 1, 1)])
	A = Int[-69 -63 -5; 119 67 1; -101 7 -88]
	@test tryindex(A, (1, 1), (2, 2), (3, 3), (4, 4), (5, 5)) == (-69, 67, -88, nothing, nothing)
	@test tryindex(A, (1, 1)) == (-69, )
	# @test tryindex(A, CartesianIndex((1, 1))) == (-69, )
	T = [(3, 1), (4, 5), (6, 2), (1, 1)]
	@test extrema_indices(T) == ((1, 6), (1, 5))
	@test extrema_indices(T...) == ((1, 6), (1, 5))
	@test append_n_times(A, 5, 1, dims = 1) == [-69 -63 -5; 119 67 1; -101 7 -88; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1]
	@test append_n_times(A, 5, 1, dims = 2) == [-69 -63 -5 1 1 1 1 1; 119 67 1 1 1 1 1 1; -101 7 -88 1 1 1 1 1]
	@test append_n_times(A, 2, 1, dims = 3) == cat([-69 -63 -5; 119 67 1; -101 7 -88], [1 1 1; 1 1 1; 1 1 1], [1 1 1; 1 1 1; 1 1 1], dims = 3)
	@test append_n_times_backwards(A, 5, 1, dims = 1) == [1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; -69 -63 -5; 119 67 1; -101 7 -88]
	@test append_n_times_backwards(A, 5, 1, dims = 2) == [1 1 1 1 1 -69 -63 -5; 1 1 1 1 1 119 67 1; 1 1 1 1 1 -101 7 -88]
	@test promote_to_nD(A, 3, 0) == cat([0 0 0; 0 0 0; 0 0 0], [-69 -63 -5; 119 67 1; -101 7 -88], [0 0 0; 0 0 0; 0 0 0], dims = 3)
	@test promote_to_3D(A, 0) == cat([0 0 0; 0 0 0; 0 0 0], [-69 -63 -5; 119 67 1; -101 7 -88], [0 0 0; 0 0 0; 0 0 0], dims = 3)
	# @test expand_as_required(A, 0, CartesianIndex((4, 4))) == [-69 -63 -5 0; 119 67 1 0; -101 7 -88 0; 0 0 0 0]
	@test expand_as_required(A, 0, [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5)]) == [-69 -63 -5 0 0; 119 67 1 0 0; -101 7 -88 0 0; 0 0 0 0 0; 0 0 0 0 0]
	# @test expand_as_required(A, 0, [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (-2, -2)]) == [-69 -63 -5 0 0; 119 67 1 0 0; -101 7 -88 0 0; 0 0 0 0 0; 0 0 0 0 0]
	@test aresamearrays(adjacencies(A, (2, 2)), [-69, 119, -101, -63, 7, -5, 1, -88])
	@test aresamearrays(adjacencies(A, (3, 3)), [67, 7, 1])
	@test aresamearrays(expanded_adjacencies(A, 0, (2, 2)), [-69, 119, -101, -63, 7, -5, 1, -88])
	@test aresamearrays(expanded_adjacencies(A, 0, (3, 3)), [67, 7, 0, 1, 0, 0, 0, 0])
	@test n_adjacent_to(A, (2, 2), 7) == 1
	@test n_adjacent_to(A, (3, 3), 7) == 1
	@test n_adjacent_to(promote_to_3D(A, 7), (3, 3, 1), 7) == 4
	@test expanded_n_adjacent_to(A, 7, (3, 3), 7) == 6
	@test expanded_n_adjacent_to(A, 7, (6, 6), 7) == 8 # (the filled element is surrounding the given index)
	@test expanded_n_adjacent_to(promote_to_3D(A, 7), 7, (3, 3, 1), 7) == 24
	@test aresamearrays(global_adjacencies_indices(A, (2, 2), 7), [(1, 1), (2, 1), (3, 1), (1, 2), (1, 3), (2, 3), (3, 3)])
	@test aresamearrays(global_adjacencies_indices(promote_to_3D(A, 7), (3, 3, 1), 7), [(2, 2, 2), (2, 3, 2), (3, 3, 2)])
	@test aresamearrays(global_adjacencies(A, (3, 3), 67), [-69, 7, 1])
	@test global_n_adjacent_to(A, (3, 3), -69, 67) == 1
    @test aresamearrays(reduce_dims(zeros(Int, 3, 3, 3)), fill(zeros(Int, 3), 3, 3))
    @test aresamearrays(reduce_dims(zeros(Int, 20, 3, 10, 100)), fill(zeros(Int, 100), 20, 3, 10))
    @test aresamearrays(reduce_dims(zeros(Int, 3, 3)), fill(zeros(Int, 3), 3))
end # end testset
