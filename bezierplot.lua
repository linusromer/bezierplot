#!/usr/bin/env lua
-- Linus Romer, published 2018 under LPPL Version 1.3c
-- version 1.3 2018-09-03
abs = math.abs
acos = math.acos
asin = math.asin
atan = math.atan
cos = math.cos
exp = math.exp
e = math.exp(1)
log = math.log
pi = math.pi
sin = math.sin
sqrt = math.sqrt
tan = math.tan
huge = math.huge

-- cube root defined for all real numbers x
function cbrt(x)
	if x < 0 then
		return -(-x)^(1/3)
	else
		return x^(1/3)
	end
end

function sgn(x)
	if x<0 then
		return -1
	elseif x>0 then
		return 1
	else
		return 0
	end
end

local function round(num, decimals)
	local result = tonumber(string.format("%." .. (decimals or 0) .. "f", num))
	if abs(result) == 0 then
		return 0
	else
		return result
	end
end

-- check if a point (x,y) satisfies xmin <= x <= xmax and ymin <= < <= ymax
local function is_in_window(x,y,xmin,xmax,ymin,ymax)
	if x >= xmin and x <= xmax and y >= ymin and y <= ymax then
		return true
	else
		return false
	end
end

local function evaluate(s)
	local tempfunc = assert(load("return " .. s))
	return tempfunc()
end


-- 5-stencil method
-- return from a graph from f in the form {{x,y},...}
-- the derivatives in form {{x,y,dy/dx,ddy/ddx,extrema,inflection},...}
local function diffgraph(func,graph,h)
	local dgraph = {}	
	local l = #graph
	if l < 4 then -- this is not worth the pain...
		for i = 1, l do 
			table.insert(dgraph,{graph[i][1],graph[i][2],0,0,0,0})
		end
	else
		local yh = func(graph[1][1]-h)
		local yhh = func(graph[1][1]-2*h)
		if yhh > -math.huge and yhh < math.huge  -- if defined at all
		and yh > -math.huge and yh < math.huge then
			dgraph[1] = {graph[1][1],graph[1][2],
			(yhh-8*yh+8*graph[2][2]-graph[3][2])/(12*h),
			(-yhh+16*yh-30*graph[1][2]+16*graph[2][2]-graph[3][2])
			/(12*h^2)}
			dgraph[2] = {graph[2][1],graph[2][2],
			(yh-8*graph[1][2]+8*graph[3][2]-graph[4][2])/(12*h),
			(-yh+16*graph[1][2]-30*graph[2][2]+16*graph[3][2]-graph[4][2])
			/(12*h^2)}
		else -- take neighbour values
			dgraph[1] = {graph[1][1],graph[1][2],
			(graph[1][2]-8*graph[2][2]+8*graph[4][2]-graph[5][2])/(12*h),
			(-graph[1][2]+16*graph[2][2]-30*graph[3][2]
			+16*graph[4][2]-graph[5][2])/(12*h^2)}
			dgraph[2] = {graph[2][1],graph[2][2],
			(graph[1][2]-8*graph[2][2]+8*graph[4][2]-graph[5][2])/(12*h),
			(-graph[1][2]+16*graph[2][2]-30*graph[3][2]
			+16*graph[4][2]-graph[5][2])/(12*h^2)}
		end
		for i = 3, l-2 do 
			table.insert(dgraph,{graph[i][1],graph[i][2],
			(graph[i-2][2]-8*graph[i-1][2]+8*graph[i+1][2]-graph[i+2][2])
			/(12*h),
			(-graph[i-2][2]+16*graph[i-1][2]-30*graph[i][2]
			+16*graph[i+1][2]-graph[i+2][2])
			/(12*h^2)})
		end
		yh = func(graph[l][1]+h)
		yhh = func(graph[l][1]+2*h)
		if yhh > -math.huge and yhh < math.huge  -- if defined at all
		and yh > -math.huge and yh < math.huge then
			dgraph[l-1] = {graph[l-1][1],graph[l-1][2],
			(graph[l-3][2]-8*graph[l-2][2]+8*graph[l][2]-yh)/(12*h),
			(-graph[l-3][2]+16*graph[l-2][2]-30*graph[l-1][2]
			+16*graph[l][2]-yh)/(12*h^2)}
			dgraph[l] = {graph[l][1],graph[l][2],
			(graph[l-2][2]-8*graph[l-1][2]+8*yh-yhh)/(12*h),
			(-graph[l-2][2]+16*graph[l-1][2]-30*graph[l][2]
			+16*yh-yhh)/(12*h^2)}
		else	
			-- take neighbour values
			dgraph[l] = {graph[l][1],graph[l][2],
			(graph[l-4][2]-8*graph[l-3][2]+8*graph[l-1][2]-graph[l][2])
			/(12*h),
			(-graph[l-4][2]+16*graph[l-3][2]-30*graph[l-2][2]
			+16*graph[l-1][2]-graph[l][2])/(12*h^2)}
			dgraph[l-1] = {graph[l-1][1],graph[l-2][2],
			(graph[l-4][2]-8*graph[l-3][2]+8*graph[l-1][2]-graph[l][2])
			/(12*h),
			(-graph[l-4][2]+16*graph[l-3][2]-30*graph[l-2][2]
			+16*graph[l-1][2]-graph[l][2])/(12*h^2)}
		end
		-- add information about being extremum / inflection point (true/false)
		for i = 1, l do 
			dgraph[i][5] = false -- dy/dx == 0 ? default, may change later
			dgraph[i][6] = false -- ddy/ddx == 0 ? default, may change later
		end
		for i = 1, l-1 do 
			-- if no gap is inbetween
			if not (dgraph[i+1][1] - dgraph[i][1] > 1.5*h) then
				-- check for dy/dx == 0 
				-- if not already determined as near dy/dx=0
				if dgraph[i][3] == 0 then
					if dgraph[i+1][3] == 0 then --take the later
						dgraph[i+1][5] = true
						dgraph[i][5] = false
					else
						dgraph[i][5] = true
					end
				elseif abs(dgraph[i][3]*dgraph[i+1][3]) 
				~= dgraph[i][3]*dgraph[i+1][3] then -- this must be near
					if abs(dgraph[i][4]) <= abs(dgraph[i+1][4]) then
						dgraph[i][5] = true
					else
						dgraph[i+1][5] = true
					end
				end
				-- check for ddy/ddx == 0 
				-- if not already determined as near ddy/ddx=0
				if (not dgraph[i][6]) and (abs(dgraph[i][4]*dgraph[i+1][4]) 
					~= dgraph[i][4]*dgraph[i+1][4]) then -- this must be near
					if abs(dgraph[i][4]) <= abs(dgraph[i+1][4]) then
						dgraph[i][6] = true
					else
						dgraph[i+1][6] = true
					end
				end
			end
		end
	end
	return dgraph
end

-- simplified diffgraph function, the function is derived only once
-- return from a graph from f in the form {{x,y},...}
-- the derivatives in form {{x,y,dy/dx},...}
-- we start with index 1 and then always jump indexjump to the next
-- index
local function diffgraphsimple(func,graph,h,indexjump)
	local dgraph = {}	
	local l = math.floor((#graph-1)/indexjump)*indexjump+1
	if l < 2 then -- this is not worth the pain...
		dgraph = {graph[1][1],graph[1][2],0}
	else
		local yh = func(graph[1][1]-h)
		local yhh = func(graph[1][1]-2*h)
		if yhh > -math.huge and yhh < math.huge  -- if defined at all
		and yh > -math.huge and yh < math.huge then
			dgraph[1] = {graph[1][1],graph[1][2],
			(yhh-8*yh+8*graph[2][2]-graph[3][2])/(12*h)}
		else -- take neighbour values
			dgraph[1] = {graph[1][1],graph[1][2],
			(graph[1][2]-8*graph[2][2]+8*graph[4][2]-graph[5][2])
			/(12*h)}
		end
		for i = 1+indexjump, l-1, indexjump do 
			table.insert(dgraph,{graph[i][1],graph[i][2],
			(graph[i-2][2]-8*graph[i-1][2]+8*graph[i+1][2]-graph[i+2][2])
			/(12*h)})
		end
		yh = func(graph[l][1]+h)
		yhh = func(graph[l][1]+2*h)
		if yhh > -math.huge and yhh < math.huge  -- if defined at all
		and yh > -math.huge and yh < math.huge then

			table.insert(dgraph,{graph[l][1],graph[l][2],
			(graph[l-2][2]-8*graph[l-1][2]+8*yh-yhh)/(12*h)})
		else	
			-- take neighbour values
			table.insert(dgraph,{graph[l][1],graph[l][2],
			(graph[l-4][2]-8*graph[l-3][2]+8*graph[l-1][2]-graph[l][2])
			/(12*h)})
		end
	end
	return dgraph
end

-- diffgraph for cubic function ax^3+bx^2+cx+d
-- returns the derivatives in form {{x,y,dy/dx,ddy/ddx},...}
-- if isinverse = true then the coordinates will be inversed
local function diffgraphcubic(graph,a,b,c,d,isinverse)
	local dgraph = {}	
	local l = #graph
	for i = 1, l do 
		if isinverse then
			dgraph[#dgraph+1] = {graph[i][2],graph[i][1],c
			+graph[i][2]*(2*b+3*a*graph[i][2]),6*a*graph[i][2]+2*b}
		else
			dgraph[#dgraph+1] = {graph[i][1],graph[i][2],c
			+graph[i][1]*(2*b+3*a*graph[i][1]),6*a*graph[i][1]+2*b}
		end
	end
	return dgraph
end

-- checks for 100 x, if the function given by funcstring
-- fits the graph g (up to maxerror) after filling in
-- the parameters a, b, c, d
-- if the graph is inverted, then isinverse has to be set true
local function do_parameters_fit(a,b,c,d,funcstring,funcgraph,maxerror,isinverse)
	if not (a > -math.huge and a < math.huge and b > -math.huge and b < math.huge  and
	c > -math.huge and c < math.huge and d > -math.huge and d < math.huge) then
		return false
	end
	local funcx = string.gsub(funcstring, "a", a)
	local funcx = string.gsub(funcx, "b", b)
	local funcx = string.gsub(funcx, "c", c)
	local funcx = string.gsub(funcx, "d", d)
	local func = assert(load("local x = ...; return "..funcx))
	for i = 1, #funcgraph, math.max(1,math.floor(0.01*#funcgraph)) do 
		if isinverse then
			if abs(func(funcgraph[i][2])-funcgraph[i][1]) 
			> maxerror then
				return false
			end
		else
			if abs(func(funcgraph[i][1])-funcgraph[i][2]) 
			> maxerror then
				return false
			end
		end
	end
	return true
end

-- f(x)=a*x^3+b*x+c
local function parameters_cubic(xp,yp,xq,yq,xr,yr,xs,ys)
	local a = (((xp^2 * xq) * yr) - ((xp^2 * xq) * ys) 
	- ((xp^2 * xr) * yq) + ((xp^2 * xr) * ys) + ((xp^2 * xs) * yq)
	- ((xp^2 * xs) * yr) - ((xp * xq^2) * yr) + ((xp * xq^2) * ys) 
	+ ((xp * xr^2) * yq) - ((xp * xr^2) * ys) - ((xp * xs^2) * yq) 
	+ ((xp * xs^2) * yr) + ((xq^2 * xr) * yp) - ((xq^2 * xr) * ys) 
	- ((xq^2 * xs) * yp) + ((xq^2 * xs) * yr) - ((xq * xr^2) * yp) 
	+ ((xq * xr^2) * ys) + ((xq * xs^2) * yp) - ((xq * xs^2) * yr) 
	+ ((xr^2 * xs) * yp) - ((xr^2 * xs) * yq) - ((xr * xs^2) * yp) 
	+ ((xr * xs^2) * yq)) / 
	(((xp^3 * xq^2) * xr) - ((xp^3 * xq^2) * xs) 
	- ((xp^3 * xq) * xr^2) + ((xp^3 * xq) * xs^2) 
	+ ((xp^3 * xr^2) * xs) - ((xp^3 * xr) * xs^2) 
	- ((xp^2 * xq^3) * xr) + ((xp^2 * xq^3) * xs) 
	+ ((xp^2 * xq) * xr^3) - ((xp^2 * xq) * xs^3) 
	- ((xp^2 * xr^3) * xs) + ((xp^2 * xr) * xs^3) 
	+ ((xp * xq^3) * xr^2) - ((xp * xq^3) * xs^2) 
	- ((xp * xq^2) * xr^3) + ((xp * xq^2) * xs^3) 
	+ ((xp * xr^3) * xs^2) - ((xp * xr^2) * xs^3) 
	- ((xq^3 * xr^2) * xs) + ((xq^3 * xr) * xs^2) 
	+ ((xq^2 * xr^3) * xs) - ((xq^2 * xr) * xs^3) 
	- ((xq * xr^3) * xs^2) + ((xq * xr^2) * xs^3))
	local b = ((((-xp^3) * xq) * yr) + ((xp^3 * xq) * ys) 
	+ ((xp^3 * xr) * yq) - ((xp^3 * xr) * ys) - ((xp^3 * xs) * yq) 
	+ ((xp^3 * xs) * yr) + ((xp * xq^3) * yr) - ((xp * xq^3) * ys) 
	- ((xp * xr^3) * yq) + ((xp * xr^3) * ys) + ((xp * xs^3) * yq) 
	- ((xp * xs^3) * yr) - ((xq^3 * xr) * yp) + ((xq^3 * xr) * ys) 
	+ ((xq^3 * xs) * yp) - ((xq^3 * xs) * yr) + ((xq * xr^3) * yp) 
	- ((xq * xr^3) * ys) - ((xq * xs^3) * yp) + ((xq * xs^3) * yr) 
	- ((xr^3 * xs) * yp) + ((xr^3 * xs) * yq) + ((xr * xs^3) * yp) 
	- ((xr * xs^3) * yq)) / 
	(((xp^3 * xq^2) * xr) - ((xp^3 * xq^2) * xs) 
	- ((xp^3 * xq) * xr^2) + ((xp^3 * xq) * xs^2) 
	+ ((xp^3 * xr^2) * xs) - ((xp^3 * xr) * xs^2) 
	- ((xp^2 * xq^3) * xr) + ((xp^2 * xq^3) * xs) 
	+ ((xp^2 * xq) * xr^3) - ((xp^2 * xq) * xs^3) 
	- ((xp^2 * xr^3) * xs) + ((xp^2 * xr) * xs^3) 
	+ ((xp * xq^3) * xr^2) - ((xp * xq^3) * xs^2) 
	- ((xp * xq^2) * xr^3) + ((xp * xq^2) * xs^3) 
	+ ((xp * xr^3) * xs^2) - ((xp * xr^2) * xs^3) 
	- ((xq^3 * xr^2) * xs) + ((xq^3 * xr) * xs^2) 
	+ ((xq^2 * xr^3) * xs) - ((xq^2 * xr) * xs^3) 
	- ((xq * xr^3) * xs^2) + ((xq * xr^2) * xs^3))
	local c = (((xp^3 * xq^2) * yr) - ((xp^3 * xq^2) * ys) 
	- ((xp^3 * xr^2) * yq) + ((xp^3 * xr^2) * ys) 
	+ ((xp^3 * xs^2) * yq) - ((xp^3 * xs^2) * yr) 
	- ((xp^2 * xq^3) * yr) + ((xp^2 * xq^3) * ys) 
	+ ((xp^2 * xr^3) * yq) - ((xp^2 * xr^3) * ys) 
	- ((xp^2 * xs^3) * yq) + ((xp^2 * xs^3) * yr) 
	+ ((xq^3 * xr^2) * yp) - ((xq^3 * xr^2) * ys) 
	- ((xq^3 * xs^2) * yp) + ((xq^3 * xs^2) * yr) 
	- ((xq^2 * xr^3) * yp) + ((xq^2 * xr^3) * ys) 
	+ ((xq^2 * xs^3) * yp) - ((xq^2 * xs^3) * yr) 
	+ ((xr^3 * xs^2) * yp) - ((xr^3 * xs^2) * yq) 
	- ((xr^2 * xs^3) * yp) + ((xr^2 * xs^3) * yq)) / 
	(((xp^3 * xq^2) * xr) - ((xp^3 * xq^2) * xs) 
	- ((xp^3 * xq) * xr^2) + ((xp^3 * xq) * xs^2) 
	+ ((xp^3 * xr^2) * xs) - ((xp^3 * xr) * xs^2) 
	- ((xp^2 * xq^3) * xr) + ((xp^2 * xq^3) * xs) 
	+ ((xp^2 * xq) * xr^3) - ((xp^2 * xq) * xs^3) 
	- ((xp^2 * xr^3) * xs) + ((xp^2 * xr) * xs^3) 
	+ ((xp * xq^3) * xr^2) - ((xp * xq^3) * xs^2) 
	- ((xp * xq^2) * xr^3) + ((xp * xq^2) * xs^3) 
	+ ((xp * xr^3) * xs^2) - ((xp * xr^2) * xs^3) 
	- ((xq^3 * xr^2) * xs) + ((xq^3 * xr) * xs^2) 
	+ ((xq^2 * xr^3) * xs) - ((xq^2 * xr) * xs^3) 
	- ((xq * xr^3) * xs^2) + ((xq * xr^2) * xs^3))
	local d = ((((xp^(3) * xq^(2)) * xr) * ys) 
	- (((xp^(3) * xq^(2)) * xs) * yr) - (((xp^(3) * xq) * xr^(2)) * ys) 
	+ (((xp^(3) * xq) * xs^(2)) * yr) + (((xp^(3) * xr^(2)) * xs) * yq) 
	- (((xp^(3) * xr) * xs^(2)) * yq) - (((xp^(2) * xq^(3)) * xr) * ys) 
	+ (((xp^(2) * xq^(3)) * xs) * yr) + (((xp^(2) * xq) * xr^(3)) * ys) 
	- (((xp^(2) * xq) * xs^(3)) * yr) - (((xp^(2) * xr^(3)) * xs) * yq) 
	+ (((xp^(2) * xr) * xs^(3)) * yq) + (((xp * xq^(3)) * xr^(2)) * ys) 
	- (((xp * xq^(3)) * xs^(2)) * yr) - (((xp * xq^(2)) * xr^(3)) * ys) 
	+ (((xp * xq^(2)) * xs^(3)) * yr) + (((xp * xr^(3)) * xs^(2)) * yq) 
	- (((xp * xr^(2)) * xs^(3)) * yq) - (((xq^(3) * xr^(2)) * xs) * yp) 
	+ (((xq^(3) * xr) * xs^(2)) * yp) + (((xq^(2) * xr^(3)) * xs) * yp) 
	- (((xq^(2) * xr) * xs^(3)) * yp) - (((xq * xr^(3)) * xs^(2)) * yp) 
	+ (((xq * xr^(2)) * xs^(3)) * yp)) / 
	(((xp^(3) * xq^(2)) * xr) - 
	((xp^(3) * xq^(2)) * xs) - ((xp^(3) * xq) * xr^(2)) 
	+ ((xp^(3) * xq) * xs^(2)) + ((xp^(3) * xr^(2)) * xs) 
	- ((xp^(3) * xr) * xs^(2)) - ((xp^(2) * xq^(3)) * xr) 
	+ ((xp^(2) * xq^(3)) * xs) + ((xp^(2) * xq) * xr^(3)) 
	- ((xp^(2) * xq) * xs^(3)) - ((xp^(2) * xr^(3)) * xs) 
	+ ((xp^(2) * xr) * xs^(3)) + ((xp * xq^(3)) * xr^(2)) 
	- ((xp * xq^(3)) * xs^(2)) - ((xp * xq^(2)) * xr^(3)) 
	+ ((xp * xq^(2)) * xs^(3)) + ((xp * xr^(3)) * xs^(2)) 
	- ((xp * xr^(2)) * xs^(3)) - ((xq^(3) * xr^(2)) * xs) 
	+ ((xq^(3) * xr) * xs^(2)) + ((xq^(2) * xr^(3)) * xs) 
	- ((xq^(2) * xr) * xs^(3)) - ((xq * xr^(3)) * xs^(2)) 
	+ ((xq * xr^(2)) * xs^(3)))
	return a, b, c, d
end

-- f(x)=a*x+b
local function parameters_affine(xp,yp,xq,yq)
	local a = (yp - yq) / (xp - xq)
	local b = ((xp * yq) - (xq * yp)) / (xp - xq)
	return a, b
end

-- what is the sum of the squared error
-- when comparing the bezier path
-- p.. control q and r .. s
-- with the graph g from index starti to endi
-- (looking at the points at roughly t=.33 and t=.67)
local function squareerror(f,g,starti,endi,qx,qy,rx,ry)
	local result = 0
	for t = .1, .9, .1 do
		x = (1-t)^3*g[starti][1]+3*t*(1-t)^2*qx+3*t^2*(1-t)*rx+t^3*g[endi][1]
		y = (1-t)^3*g[starti][2]+3*t*(1-t)^2*qy+3*t^2*(1-t)*ry+t^3*g[endi][2]
		result = result + (y-f(x))^2
	end
	return result
end

-- converts a table with bezier point information
-- to a string with rounded values 
-- the path is reversed, if rev is true
-- e.g. if bezierpoints = {{0,1},{2,3,4,5,6,7},{8,9,10,11,12,13}}
-- then 
-- (0,1) .. controls (2,3) and (4,5) .. (6,7) .. controls 
-- (8,9) and (10,11) .. (12,13)
-- will be returned
-- the notation "pgfplots" will change the notation to
-- YES: \addplot coordinates {(0,1) (6,7) (2,3) (4,5) (6,7) (12,13) (8,9) (10,11)}
-- NO: 0  1 \\ 6  7 \\ 2  3 \\ 4  5 \\ \\ 6  7 \\ 12 13 \\ 8  9 \\ 10 11 \\
-- As pgfplots does not connect the bezier segments
-- reverse paths are not implemented 
local function beziertabletostring(bezierpoints,rndx,rndy,rev,notation)
	local bezierstring = ""
	local b = {{round(bezierpoints[1][1],rndx),round(bezierpoints[1][2],rndy)}} -- rounded and then 
	-- reduced points (if identical after rounding)
	-- rounding
	for i = 2, #bezierpoints do
		-- check if x--coordinates are identical
		if round(bezierpoints[i][#bezierpoints[i]-1],rndx) ~=  b[#b][#b[#b]-1] then
			b[#b+1] = {}
			for j = 1, #bezierpoints[i] do
				if j % 2 == 0 then -- x coordinate
					b[#b][j] = round(bezierpoints[i][j],rndx)
				else
					b[#b][j] = round(bezierpoints[i][j],rndy)
				end
			end
		end
	end
	if #b > 1 then -- if not empty or single point
		-- check if bezierstring contains only straight lines
		local onlystraightlines = true
		for i = 1, #b do
			if #b[i] > 2 then
				onlystraightlines = false
				break
			end
		end
		if onlystraightlines then
			if rev then
				bezierstring = "(" .. b[#b][1] .. "," .. b[#b][2] ..")"
					for i = #b-1, 1, -1 do
						bezierstring = bezierstring .. 
							" -- (" .. b[i][1] .. "," .. b[i][2] ..")"
					end
			else
				if notation == "pgfplots" then
					bezierstring = "\\addplot coordinates {(" 
						.. b[1][1] .. "," .. b[1][2] .. ") (" 
						.. b[2][1] .. "," .. b[2][2] .. ") (" 
						.. b[1][1] .. "," .. b[1][2] .. ") (" 
						.. b[2][1] .. "," .. b[2][2] .. ") }" 
				else -- notation = tikz
					bezierstring = "(" .. b[1][1] .. "," .. b[1][2] ..")"
					for i = 2, #b do
						bezierstring = bezierstring .. 
							" -- (" .. b[i][1] .. "," .. b[i][2] ..")"
					end
				end	
			end
		else
			if rev then
				bezierstring = "(" .. b[#b][#b[#b]-1] .. "," 
				.. b[#b][#b[#b]] ..")" -- initial point
				for i = #b, 2, -1 do
					if #b[i] >= 6 then -- cubic bezier spline
						bezierstring = bezierstring .. " .. controls (" 
						.. b[i][3] .. "," .. b[i][4] ..") and ("
						.. b[i][1] .. "," .. b[i][2] .. ") .. (" 
						.. b[i-1][#b[i-1]-1] .. "," .. b[i-1][#b[i-1]]..")"
					else
						bezierstring = bezierstring .. " (" 
						.. b[i-1][#b[i-1]-1] .. "," .. b[i-1][#b[i-1]] ..")"
					end
				end
			else
				if notation == "pgfplots" then
					bezierstring = "\\addplot coordinates {"
					for i = 1, #b-1 do
						if #b[i+1] >= 6 then -- cubic bezier spline
							bezierstring = bezierstring .. "("
							.. b[i][#b[i]-1] .. "," .. b[i][#b[i]] .. ") (" 
							.. b[i+1][5] .. "," .. b[i+1][6] .. ") ("  
							.. b[i+1][1] .. "," .. b[i+1][2] .. ") (" 
							.. b[i+1][3] .. "," .. b[i+1][4] .. ") " 
						end
					end
					bezierstring = bezierstring .. "}"
				else -- notation = tikz
					bezierstring = "(" .. b[1][1] .. "," 
					.. b[1][2] ..")" -- initial point
					for i = 2, #b do
						if #b[i] >= 6 then -- cubic bezier spline
							bezierstring = bezierstring .. " .. controls (" 
							.. b[i][1] .. "," .. b[i][2] ..") and ("
							.. b[i][3] .. "," .. b[i][4] .. ") .. (" 
							.. b[i][5] .. "," .. b[i][6]..")"
						else
							bezierstring = bezierstring .. " (" 
							.. b[i][1] .. "," .. b[i][2] ..")"
						end
					end
				end
			end
		end
	end
	return bezierstring
end

-- take end points of a graph g of the function f
-- (from indices starti to endi)
-- without extrema or inflection points inbetween 
-- and try to approximate it with a cubic bezier curve
-- (round to rndx and rndy when printing)
-- if maxerror <= 0, the function will not be recursive anymore
local function graphtobezierapprox(f,g,starti,endi,maxerror)
	local px = g[starti][1]
	local py = g[starti][2]
	local dp = g[starti][3]
	local sx = g[endi][1]
	local sy = g[endi][2]
	local ds = g[endi][3]
	-- we compute the corner point c, where the controls would meet
	local cx = ((dp * px) - (ds * sx) - py + sy) / (dp - ds)
	local cy = (dp * ((ds * px) - (ds * sx) - py + sy) / (dp - ds)) + py
	-- now we slide q between p and c & r between s and c
	-- and search for the best qx and best rx
	local qx = px+.01*(cx-px)
	local qy = py+.01*(cy-py)
	local rx = sx+.01*(cx-sx)
	local ry = sy+.01*(cy-sy)
	local err = squareerror(f,g,starti,endi,qx,qy,rx,ry)
	for i = 2, 99 do
		for j = 2, 99 do
			xa = px+i*.01*(cx-px)
			ya = py+i*.01*(cy-py)
			xb = sx+j*.01*(cx-sx)
			yb = sy+j*.01*(cy-sy)
			-- now check, if xa and xb fit better
			-- than the last qx and rx did
			-- (sum of squares must be smaller)
			local newerror = squareerror(f,g,starti,endi,xa,ya,xb,yb)
			if newerror < err then
				qx = xa
				qy = ya
				rx = xb
				ry = yb
				err = newerror
			end
		end
	end
	if maxerror > 0 then
		-- check if it is close enough: (recycling err, xa, ya)
		err = 0
		for t = .1, .9, .1 do
			xa = (1-t)^3*g[starti][1]+3*t*(1-t)^2*qx+3*t^2*(1-t)*rx+t^3*g[endi][1]
			ya = (1-t)^3*g[starti][2]+3*t*(1-t)^2*qy+3*t^2*(1-t)*ry+t^3*g[endi][2]
			if abs(ya-f(xa)) > err then
				err = abs(ya-f(xa))
				err = abs(ya-f(xa))
			end
		end
		if err <= maxerror then
			return {qx,qy,rx,ry,sx,sy}
		else
			-- search for an intermediate point where the graph has the same
			-- slope as the line from the start point to the end point:
			local interindex = math.floor(.5*starti+.5*endi) -- will change
			for i = starti + 1, endi - 1 do
				if abs(g[i][3]-(g[endi][2]-g[starti][2])
				/(g[endi][1]-g[starti][1])) 
				< abs(g[interindex][3]-(g[endi][2]-g[starti][2])
				/(g[endi][1]-g[starti][1])) then
					interindex = i
				end
			end
			local left = graphtobezierapprox(f,g,starti,interindex,maxerror)
			local right = graphtobezierapprox(f,g,interindex,endi,maxerror)
			for i=1, #right do --now append the right to the left:
				left[#left+1] = right[i]
			end
			return left
		end
	else
		return {qx,qy,rx,ry,sx,sy}
	end
end

-- like above but exact for quadratic and cubic (if not inverse)
-- resp. exact for squareroot and cuberoot (if inverse)
local function graphtobezier(g,starti,endi,isinverse)
	local px = g[starti][1]
	local py = g[starti][2]
	local dp = g[starti][3]
	local sx = g[endi][1]
	local sy = g[endi][2]
	local ds = g[endi][3]
	local qx = px+(sx-px)/3
	local rx = px+2*(sx-px)/3
	local qy = py+(qx-px)*dp
	local ry = sy+(rx-sx)*ds
	if isinverse then
		return {qy,qx,ry,rx,sy,sx}
	else
		return {qx,qy,rx,ry,sx,sy}
	end
end

-- just for debugging:
local function printtable(t)
	for i = 1,#t do
		for j = 1, #t[i] do
			io.write(t[i][j].." ")
		end
		io.write("\n")
	end
end

-- main function
function bezierplot(functionstring,xminstring,xmaxstring,yminstring,ymaxstring,samplesstring,notation)
	local fstringreplaced = string.gsub(functionstring, "%*%*", "^")
	local f = assert(load("local x = ...; return " .. fstringreplaced)) 
	local xmin = evaluate(xminstring)
	local xmax = evaluate(xmaxstring)
	local ymin = evaluate(yminstring)
	local ymax = evaluate(ymaxstring)
	local samples = evaluate(samplesstring)
	local isreverse = false
	if xmin > xmax then
		isreverse = true
	elseif xmin == xmax then
		xmax = xmin + 10
	end
	xmin, xmax = math.min(xmin,xmax), math.max(xmin,xmax)
	if ymin == ymax then
		ymax = ymin + 10
	end
	ymin, ymax = math.min(ymin,ymax), math.max(ymin,ymax)
	local xsteps = 50000 
	-- if samples < 2 the samples will be chosen as wisely as possible
	local arbitrary_samples = true
	if samples >= 2 then
		arbitrary_samples = false
		xsteps = (samples-1)*math.max(2,math.floor(xsteps/(samples-1)))
	end
	local xstep = (xmax-xmin)/xsteps
	-- the output of the x coordinates will be rounded to rndx digits
	local rndx = math.max(0,math.floor(5.5-log(xmax-xmin)/log(10)))
	local xerror = abs(xmax-xmin)/(10^rndx)
	-- the output of the y coordinates will be rounded to rndy digits
	local rndy = math.max(0,math.floor(5.5-log(ymax-ymin)/log(10)))
	local yerror = (ymax-ymin)/(10^rndy)
	-- determine parts of the graph that are inside window
	local graphs = {} -- graph split to the connected parts
	local graph = {} -- graphs concatenated (needed for function type)
	local outside = true -- value is outside window
	local i = 0
	local j = 0
	local yminreal -- determine the real minimimum of the y coord.
	local ymaxreal -- just decring
	local yminrealfound = false
	local ymaxrealfound = false
	for n = 0, xsteps do
		local x = xmin + n/xsteps*(xmax-xmin)
		if n == xsteps then
			x = xmax
		end
		local y = f(x)
		if (y >= ymin-.1*yerror and ymin ~= -huge or y > ymin and ymin == -huge)
		and (y <= ymax+.1*yerror and ymax ~= huge or y < ymax and ymax == huge)
		then -- inside
			if outside then -- if it was outside before
				outside = false
				j = 0
				i = i + 1
				graphs[i] = {}
			end
			j = j + 1
			graphs[i][j] = {x,y}
			graph[#graph+1] = {x,y}
			if not yminrealfound or yminrealfound and y < yminreal then
				yminreal = y
				yminrealfound = true
			end
			if not ymaxrealfound or ymaxrealfound and y > ymaxreal then
				ymaxreal = y
				ymaxrealfound = true
			end
		else
			outside = true
		end
	end
	
	-- some redefinitions
	if #graph ~= 0 and yminreal ~= ymaxreal then
		ymin = yminreal
		ymax = ymaxreal
		rndy = math.max(0,math.floor(5.5-log(ymax-ymin)/log(10)))
		yerror = (ymax-ymin)/(10^rndy)
	end
	
	-- check for the function type (for this, we need the concatenated
	-- parts of the graph)
	-- go through the connected parts
	local functiontype = "unknown"
	local a, b, c, d -- possible function parameter
	-- check for affine functions:
	local l = #graph
	a, b = parameters_affine(graph[1][1],graph[1][2],
	graph[l][1],graph[l][2])
	if do_parameters_fit(a,b,0,0,"a*x+b",graph,yerror,false) then
		functiontype = "affine"
	else -- check for cubic functions (includes quadratic functions)
		a, b, c, d = parameters_cubic(graph[1][1],graph[1][2],
		graph[math.floor(l/3)][1],graph[math.floor(l/3)][2],
		graph[math.floor(2*l/3)][1],graph[math.floor(2*l/3)][2],
		graph[l][1],graph[l][2])
		if do_parameters_fit(a,b,c,d,"a*x^3+b*x^2+c*x+d",graph,
		yerror,false) then
			functiontype = "cubic"
		else -- check for cuberoot functions (includes squareroots)
			a, b, c, d = parameters_cubic(graph[1][2],graph[1][1],
			graph[math.floor(l/3)][2],graph[math.floor(l/3)][1],
			graph[math.floor(2*l/3)][2],graph[math.floor(2*l/3)][1],
			graph[l][2],graph[l][1])
			if do_parameters_fit(a,b,c,d,"a*x^3+b*x^2+c*x+d",graph,
			xerror,true) then
				functiontype = "cuberoot"
			end
		end
	end
			
	local bezierpoints = {}
	-- the bezier path (0,1) .. controls 
	-- (2,3) and (4,5) .. (6,7) .. controls 
	-- (8,9) and (10,11) .. (12,13)
	-- will be stored as
	-- bezierpoints={{0,1},{2,3,4,5,6,7},{8,9,10,11,12,13}}
	
	if functiontype == "affine" then
		if arbitrary_samples then
			bezierpoints = {{graph[1][1],graph[1][2]},{graph[#graph][1],
			graph[#graph][2]}}
		else -- we can here savely assume that graphs has only one part,
		-- therefore graphs[1]=graph
			for i = 1, #graph, math.floor(xsteps/(samples-1)) do
				bezierpoints[#bezierpoints+1] = {graph[i][1],graph[i][2]}
			end
		end
	elseif functiontype == "cubic" then 
		local extrema_inflections = {} -- store the extrema and
		-- inflection points for arbitrary samples
		if arbitrary_samples then
			if math.abs(a) < yerror*1e-10 then -- quadratic case (one extremum)
				if is_in_window(-c/(2*b),(-c^2+4*b*d)/(4*b),xmin,xmax,
				ymin,ymax) then
					extrema_inflections = {{-c/(2*b),(-c^2+4*b*d)/(4*b)}}
				end
			else -- cubic case (two extrema and one inflection point)
				-- we order the points with the help of sgn
				-- check for first extrema
				if is_in_window((-sgn(a)*sqrt(-3*a*c+b^2)-b)/(3*a),
				(2*b^3+27*a^2*d-9*a*b*c+sqrt(b^2-3*a*c)*sgn(a)*
				(2*b^2-6*a*c))/(27*a^2),xmin,xmax,ymin,ymax) then
					extrema_inflections[#extrema_inflections+1] = 
					{(-sgn(a)*sqrt(-3*a*c+b^2)-b)/(3*a),(2*b^3+27*a^2*d-
					9*a*b*c+sqrt(b^2-3*a*c)*sgn(a)*(2*b^2-6*a*c))/(27*a^2)}
				end
				-- check for inflection point (has to be inbetween)
				if is_in_window(-b/(3*a),(2*b^3+27*a^2*d-9*a*b*c)
				/(27*a^2),xmin,xmax,ymin,ymax) then
					extrema_inflections[#extrema_inflections+1]={-b/(3*a),
					(2*b^3+27*a^2*d-9*a*b*c)/(27*a^2)}
				end
				-- check for second extrema
				if is_in_window((sgn(a)*sqrt(-3*a*c+b^2)-b)/(3*a),
				(2*b^3+27*a^2*d-9*a*b*c+sqrt(b^2-3*a*c)*sgn(a)*
				(-2*b^2+6*a*c))/(27*a^2),xmin,xmax,ymin,ymax) then
					extrema_inflections[#extrema_inflections+1] = 
					{(sgn(a)*sqrt(-3*a*c+b^2)-b)/(3*a),(2*b^3+27*a^2*d-
					9*a*b*c+sqrt(b^2-3*a*c)*sgn(a)*(-2*b^2+6*a*c))/(27*a^2)}
				end
			end
		end
		for part = 1, #graphs do 
			bezierpoints[#bezierpoints+1] = {graphs[part][1][1],
			graphs[part][1][2]} -- initial points
			local graphsamples = {}-- will be the graph reduced to the 
			-- samples (or the most important points)
			local dg -- will be the differentiated graph
			if arbitrary_samples then -- add extrema and inflection 
			-- points to the border points
				graphsamples = {{graphs[part][1][1],
					graphs[part][1][2]}}
				for j = 1, #extrema_inflections do
					if extrema_inflections[j][1] > math.min(
					graphs[part][1][1] ,graphs[part][#graphs[part]][1])
					and extrema_inflections[j][1] < math.max(
					graphs[part][1][1] ,graphs[part][#graphs[part]][1])
					then
						graphsamples[#graphsamples+1] = 
							{extrema_inflections[j][1],
							extrema_inflections[j][2]}
					end
				end
				graphsamples[#graphsamples+1] = 
					{graphs[part][#graphs[part]][1],
					graphs[part][#graphs[part]][2]}
			else
				for i = 1, #graphs[part], xsteps/(samples-1) do
					graphsamples[#graphsamples+1] = 
						{graphs[part][i][1],graphs[part][i][2]}
				end
			end
			dg = diffgraphcubic(graphsamples,a,b,c,d,false)
			for i = 2, #dg do
				bezierpoints[#bezierpoints+1] = graphtobezier(dg,i-1,i,false)
			end
		end
	elseif functiontype == "cuberoot" then 
		local inflection = {} -- store the inflection point
		if arbitrary_samples and math.abs(a) ~= 0
		and is_in_window((2*b^3+27*a^2*d-9*a*b*c)/(27*a^2),-b/(3*a),
		xmin,xmax,ymin,ymax) then
			inflection = {(2*b^3+27*a^2*d-9*a*b*c)/(27*a^2),-b/(3*a)}
		end
		-- (there cannot be more than one part)
		bezierpoints[#bezierpoints+1] = {graphs[1][1][1],
		graphs[1][1][2]} -- initial points
		local graphsamples = {}-- will be the graph reduced to the 
		-- samples (or the most important points)
		local dg -- will be the differentiated graph
		if arbitrary_samples then -- add inflection point (if exis.)
			graphsamples = {{graphs[1][1][1],
				graphs[1][1][2]}}
			if #inflection > 0 and inflection[1] > math.min(
				graphs[1][1][1],graphs[1][#graphs[1]][1])
				and inflection[1] < math.max(
				graphs[1][1][1],graphs[1][#graphs[1]][1])
				then
					graphsamples[#graphsamples+1] = 
						{inflection[1],inflection[2]}
			end
			graphsamples[#graphsamples+1] = 
				{graphs[1][#graphs[1]][1],
				graphs[1][#graphs[1]][2]}
		else
			for i = 1, #graphs[1], xsteps/(samples-1) do
				graphsamples[#graphsamples+1] = 
					{graphs[1][i][1],graphs[1][i][2]}
			end
		end
		dg = diffgraphcubic(graphsamples,a,b,c,d,true)
		for i = 2, #dg do
			bezierpoints[#bezierpoints+1] = graphtobezier(dg,i-1,i,true)
		end
	else	
	---------- generic case (no special function) ----------------		
		if arbitrary_samples then
			-- go through the connected parts
			for part = 1, #graphs do 
				local dg = diffgraph(f,graphs[part],xstep)
				bezierpoints[#bezierpoints+1] = {dg[1][1],dg[1][2]}
				local startindex = 1
				for k = 2, #dg do
					if dg[k][5] or dg[k][6] then -- extrema and inflection points
						local tobeadded = graphtobezierapprox(
						f,dg,startindex,k,10*yerror)
						-- tobeadded may contain a multiple of 6 entries
					-- e.g. {1,2,3,4,5,6,7,8,9,10,11,12}
						for i = 1, math.floor(#tobeadded/6) do
							bezierpoints[#bezierpoints+1] = {}
							for j = 1, 6 do
								bezierpoints[#bezierpoints][j] = tobeadded[(i-1)*6+j]
							end
						end
						startindex = k
					end
				end
				if startindex ~= #dg then -- if no special points inbetween
					local tobeadded = graphtobezierapprox(f,dg,
					startindex,#dg,10*yerror)
					-- tobeadded may contain a multiple of 6 entries
					-- e.g. {1,2,3,4,5,6,7,8,9,10,11,12}
					for i = 1, math.floor(#tobeadded/6) do
						bezierpoints[#bezierpoints+1] = {}
						for j = 1, 6 do
							bezierpoints[#bezierpoints][j] = tobeadded[(i-1)*6+j]
						end
					end
				end
			end
		else -- fixed samples in the generic case
			-- go through the connected parts
			for part = 1, #graphs do 
				local dg = diffgraphsimple(f,graphs[part],xstep,
					math.floor(0.5+xsteps/(samples-1)))
				bezierpoints[#bezierpoints+1] = {dg[1][1],dg[1][2]} -- initial points
				for i = 2, #dg do
					bezierpoints[#bezierpoints+1] = graphtobezier(dg,i-1,i,false)
				end
			end
		end
	end
	return beziertabletostring(bezierpoints,rndx,rndy,isreverse,notation)		
end

-- main program --

if not pcall(debug.getlocal, 4, 1) then
--if debug.getinfo(3) == nil then
	if #arg >= 1 then
		local xmin = -5
		local xmax = 5
		if #arg >= 2 then 
			xmin = arg[2]
		end
		if #arg >= 3 then
			xmax = arg[3]
		end
		local ymin = -5
		local ymax = 5
		if #arg >= 4 then 
			ymin = arg[4]
		end
		if #arg >= 5 then 
			ymax = arg[5]
		end
		local samples = 0
		if #arg >= 6 then 
			samples = arg[6]
		end
		local notation = "tikz"
		if #arg >= 7 then 
			notation = arg[7] 
		end
		print(bezierplot(arg[1],xmin,xmax,ymin,ymax,samples,notation))
	end
end



