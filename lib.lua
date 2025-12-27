local STRING_PADDING = 15

function center_pad(str, total_width, pad_char)
	pad_char = pad_char or " "

	local str_len = #str

	if str_len >= total_width then
		return str
	end

	local total_pad = total_width - str_len

	local left_pad = math.floor(total_pad / 2)
	local right_pad = total_pad - left_pad

	return string.rep(pad_char, left_pad) .. str .. string.rep(pad_char, right_pad)
end

local linearAlgebra = {}
linearAlgebra.Matrix = {}

local abs = math.abs
local floor = math.floor
local sqrt = math.sqrt

local cos = math.cos
local sin = math.sin

local table_create = table.create
local table_insert = table.insert

--[[ ROW CLASS ]]
--
local _rowClass = {}
local _newRow = function(values)
	local instance = setmetatable({}, _rowClass)
	rawset(instance, "_values", values)
	return instance
end

_rowClass.__index = function(row, key)
	if type(key) == "number" then
		return row._values[key]
	elseif key == "Length" then
		return #row._values
	end
	error("Unrecognized key: " .. tostring(key))
end

_rowClass.__newindex = function(row, key, value)
	if type(key) == "number" then
		if key < 1 or key > #row._values then
			error("Row index out of bounds", 2)
		end
		row._values[key] = value
	else
		error("Cannot set property '" .. tostring(key) .. "' on row", 2)
	end
end

_rowClass.__tostring = function(row)
	local result = table_create(row.Length)
	for i = 1, row.Length do
		result[i] = center_pad(`[{i}]:{tostring(row[i])}`, STRING_PADDING)
	end
	return table.concat(result, " ")
end

_rowClass.__mul = function(left, right)
	if type(left) == "number" or type(right) == "number" then
		local row = type(left) == "number" and right or left
		local scalar = type(left) == "number" and left or right

		if row.Length == 1 then
			local leftScalar = type(left) == "number" and left or left[1]
			local rightScalar = type(right) == "number" and right or right[1]
			return leftScalar * rightScalar
		else
			local resultArray = {}

			for i = 1, row.Length do
				resultArray[i] = scalar * row[i]
			end

			return _newRow(resultArray)
		end
	elseif left.Length == 1 or right.Length == 1 then
		if left.Length == 1 and right.Length == 1 then
			return left[1] * right[1]
		elseif left.Length == 1 then
			return left[1] * right
		else
			return left * right[1]
		end
	end

	error("Multiplication of two rows is undefined", 2)
end

local _matrixClass = {}
local _vectorClass = {}

local _newMatrix = function(rows)
	local rowInstances = table_create(#rows)
	for i = 1, #rows do
		rowInstances[i] = _newRow(rows[i])
	end

	local instance = setmetatable({}, _matrixClass)
	rawset(instance, "_rows", rowInstances)
	return instance
end

_matrixClass.__index = function(mat, key)
	if type(key) == "number" then
		return mat._rows[key]
	elseif key == "Shape" then
		return { #mat._rows, mat._rows[1].Length }
	elseif _matrixClass[key] then
		return _matrixClass[key]
	end
	error("Unrecognized key: " .. tostring(key))
end

_matrixClass.__newindex = function(mat, key, value)
	error(
		"Cannot assign to matrix property '" .. tostring(key) .. "'. Use mat[row][col] = value to modify elements.",
		2
	)
end

_matrixClass.__tostring = function(mat)
	local result = table_create(mat.Shape[1])
	for i = 1, mat.Shape[1] do
		result[i] = center_pad(`[{i}]: {tostring(mat[i])}`, STRING_PADDING)
	end
	return table.concat(result, "\n")
end

_matrixClass.__mul = function(left, right)
	local resultRows = {}

	if type(left) == "number" or type(right) == "number" then
		local mat = type(left) == "number" and right or left
		local scalar = type(left) == "number" and left or right

		for i = 1, mat.Shape[1] do
			resultRows[i] = {}
			for j = 1, mat.Shape[2] do
				resultRows[i][j] = scalar * mat[i][j]
			end
		end
	else
		if left.Shape[2] ~= right.Shape[1] then
			error(
				"Cannot multiply matrices of sizes ("
					.. left.Shape[1]
					.. " x "
					.. left.Shape[2]
					.. ") and ("
					.. right.Shape[1]
					.. " x "
					.. right.Shape[2]
					.. ")",
				2
			)
		end

		for i = 1, left.Shape[1] do
			resultRows[i] = {}
			for j = 1, right.Shape[2] do
				resultRows[i][j] = 0
				for k = 1, left.Shape[2] do
					resultRows[i][j] = resultRows[i][j] + (left[i][k] * right[k][j])
				end
			end
		end
	end

	return _newMatrix(resultRows)
end

function _matrixClass:add(i, j, v)
	self._rows[i + 1]._values[j + 1] += v

	return self
end

linearAlgebra.Matrix.new = function(rows)
	return _newMatrix(rows)
end

linearAlgebra.Matrix.identity = function(n)
	if n <= 0 then
		error("Matrix dimension must be positive", 2)
	end

	local rows = table_create(n)
	for i = 1, n do
		local row = table_create(n)
		for j = 1, n do
			row[j] = (i == j) and 1 or 0
		end
		rows[i] = row
	end
	return _newMatrix(rows)
end

linearAlgebra.Matrix.zeros = function(m, n)
	if m <= 0 or n <= 0 then
		error("Matrix dimensions must be positive", 2)
	end

	local rows = table_create(m)
	for i = 1, m do
		local row = table_create(n)
		for j = 1, n do
			row[j] = 0
		end
		rows[i] = row
	end
	return _newMatrix(rows)
end

linearAlgebra.Matrix.diagonal = function(values)
	local n = #values
	if n == 0 then
		error("Cannot create diagonal matrix from empty array", 2)
	end

	local resultRows = table_create(n)
	for i = 1, n do
		local row = table_create(n)
		for j = 1, n do
			row[j] = (i == j) and values[i] or 0
		end
		resultRows[i] = row
	end
	return _newMatrix(resultRows)
end

function _matrixClass:isDiagonal()
	local m, n = self.Shape[1], self.Shape[2]
	if m ~= n then
		return false
	end

	for i = 1, m do
		for j = 1, n do
			if i ~= j and self[i][j] ~= 0 then
				return false
			end
		end
	end
	return true
end

function _matrixClass:lu()
	if self.Shape[1] ~= self.Shape[2] then
		error("LU decomposition requires a square matrix", 2)
	end

	local n = self.Shape[1]
	if n == 0 then
		error("Cannot decompose empty matrix", 2)
	end

	local L = table_create(n)
	local U = table_create(n)
	local P = table_create(n)
	local sign = 1

	for i = 1, n do
		L[i] = table_create(n)
		U[i] = table_create(n)
		P[i] = table_create(n)
		for j = 1, n do
			L[i][j] = (i == j) and 1 or 0
			U[i][j] = self[i][j]
			P[i][j] = (i == j) and 1 or 0
		end
	end

	for k = 1, n do
		local maxVal = abs(U[k][k])
		local maxRow = k
		for i = k + 1, n do
			local val = abs(U[i][k])
			if val > maxVal then
				maxVal = val
				maxRow = i
			end
		end

		if maxRow ~= k then
			U[k], U[maxRow] = U[maxRow], U[k]
			P[k], P[maxRow] = P[maxRow], P[k]
			if k > 1 then
				for j = 1, k - 1 do
					L[k][j], L[maxRow][j] = L[maxRow][j], L[k][j]
				end
			end
			sign = -sign
		end

		if abs(U[k][k]) < 1e-12 then
			error("Matrix is singular or nearly singular", 2)
		end

		local pivot = U[k][k]
		for i = k + 1, n do
			local factor = U[i][k] / pivot
			L[i][k] = factor
			U[i][k] = 0
			for j = k + 1, n do
				U[i][j] = U[i][j] - factor * U[k][j]
			end
		end
	end

	return {
		L = _newMatrix(L),
		U = _newMatrix(U),
		P = _newMatrix(P),
		sign = sign,
	}
end

function _matrixClass:solve(b)
	if self.Shape[1] ~= self.Shape[2] then
		error("Coefficient matrix must be square", 2)
	end
	if b.Shape[2] ~= 1 then
		error("Right-hand side must be a column vector", 2)
	end
	if self.Shape[1] ~= b.Shape[1] then
		error("Matrix and vector dimensions do not match", 2)
	end

	local n = self.Shape[1]
	local decomp = self:lu()
	local L, U, P = decomp.L, decomp.U, decomp.P

	local Pb = P * setmetatable(b, _matrixClass)

	local y = table_create(n)
	for i = 1, n do
		local sum = 0
		for j = 1, i - 1 do
			sum = sum + L[i][j] * y[j]
		end
		y[i] = Pb[i][1] - sum
	end

	local x = table_create(n)
	for i = n, 1, -1 do
		local sum = 0
		for j = i + 1, n do
			sum = sum + U[i][j] * x[j]
		end
		x[i] = (y[i] - sum) / U[i][i]
	end

	local xRows = table_create(n)
	for i = 1, n do
		xRows[i] = { x[i] }
	end
	return setmetatable(_newMatrix(xRows), _vectorClass)
end

function _matrixClass:det()
	if self.Shape[1] ~= self.Shape[2] then
		error("Determinant requires a square matrix", 2)
	end

	local decomp = self:lu()
	local U = decomp.U
	local sign = decomp.sign
	local n = self.Shape[1]

	local det = 1
	for i = 1, n do
		det = det * U[i][i]
	end

	return sign * det
end

function _matrixClass:set(i, j, value)
	local m, n = self.Shape[1], self.Shape[2]
	if i < 1 or i > m or j < 1 or j > n then
		error(string.format("Index out of bounds: (%d, %d) for a %d x %d matrix", i, j, m, n), 2)
	end

	self[i + 1][j + 1] = value
	return self
end

-- VECTOR CLASS

_vectorClass.__index = function(vec, key)
	if type(key) == "number" then
		return vec._rows[key][1]
	elseif key == "Shape" then
		return { #vec._rows, 1 }
	elseif _vectorClass[key] then
		return _vectorClass[key]
	end
	error("Unrecognized key: " .. tostring(key))
end

_vectorClass.__newindex = function(self, key, val)
	self._rows[key][1] = val
end

_vectorClass.__tostring = function(vec)
	local result = table_create(vec.Shape[1])
	for i = 1, vec.Shape[1] do
		result[i] = `[{i}]: {tostring(vec[i])}`
	end
	return table.concat(result, "\n")
end

function _vectorClass:add(i, v)
	self[i + 1] += v

	return self
end

function _vectorClass:set(i, value)
	if self.Shape[2] ~= 1 then
		error("Vector.set can only be used on column vectors", 2)
	end
	self[i + 1] = value
	return self
end

linearAlgebra.Vector = {}

linearAlgebra.Vector.new = function(values)
	local n = #values
	local rows = table_create(n)
	for i = 1, n do
		rows[i] = { values[i] }
	end
	return setmetatable(_newMatrix(rows), _vectorClass)
end

linearAlgebra.Vector.zeros = function(n)
	return setmetatable(linearAlgebra.Matrix.zeros(n, 1), _vectorClass)
end

return linearAlgebra
