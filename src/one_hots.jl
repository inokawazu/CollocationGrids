function one_hots(n, order=1)
	order == 1 && return [[Int(i==j) for j in 1:n] for i in 1:n]
	vec = [one_hots(n, 1) for _ in 1:order]
	map(sum, Iterators.product(vec...))
end
