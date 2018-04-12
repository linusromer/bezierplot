#!/usr/bin/env lua
-- Linus Romer, published 2018 under LPPL Version 1.3c
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

function round(num, decimals)
	local result = tonumber(string.format("%." .. (decimals or 0) .. "f", num))
	if abs(result) == 0 then
		return 0
	else
		return result
	end
end

-- 5-stencil method
-- return from a graph from f in the form {{x,y},...}
-- the derivatives in form {{x,y,dy/dx,ddy/ddx},...}
function diffgraph(func,graph,h)
	local dgraph = {}	
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
	local l = #graph
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
			if not dgraph[i][5] then 
				if dgraph[i][3] == 0 then
					dgraph[i][5] = true
				elseif abs(dgraph[i][3]*dgraph[i+1][3]) 
				~= dgraph[i][3]*dgraph[i+1][3] then -- this must be near
					if abs(dgraph[i][4]) <= abs(dgraph[i+1][4]) then
						dgraph[i][5] = true
					else
						dgraph[i+1][5] = true
					end
				end
			end
			-- check for ddy/ddx == 0 
			-- if not already determined as near ddy/ddx=0
			if not dgraph[i][6] then 
				if abs(dgraph[i][4]*dgraph[i+1][4]) 
				~= dgraph[i][4]*dgraph[i+1][4] then -- this must be near
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

-- checks for 100 x, if the function given by funcstring
-- fits the graph g (up to maxerror) after filling in
-- the parameters a, b, c, d
-- if the graph is inverted, then isinverse has to be set true
function do_parameters_fit(a,b,c,d,funcstring,funcgraph,maxerror,isinverse)
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
function parameters_cubic(xp,yp,xq,yq,xr,yr,xs,ys)
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
function parameters_affine(xp,yp,xq,yq)
	local a = (yp - yq) / (xp - xq)
	local b = ((xp * yq) - (xq * yp)) / (xp - xq)
	return a, b
end

-- returns true iff the function is of type f(x)=a*x^3+b*x^2+c*x+d
-- a, b, c, d being real numbers
function is_cubic(graph,maxerror)
	local l = #graph
	local a, b, c, d = parameters_cubic(graph[1][1],graph[1][2],
	graph[math.floor(l/3)][1],graph[math.floor(l/3)][2],
	graph[math.floor(2*l/3)][1],graph[math.floor(2*l/3)][2],
	graph[l][1],graph[l][2])
	return do_parameters_fit(a,b,c,d,"a*x^3+b*x^2+c*x+d",graph,
	maxerror,false)
end

-- returns true iff the function is of type f(x)=a*x^3+b*x^2+c*x+d
-- a, b, c, d being real numbers
-- this takes several graph parts
-- the idea is to have a possibility to avoid tan(x)
function are_cubic(graphs,maxerror)
	if is_cubic(graphs[1],maxerror) then
		if #graphs < 2 then
			return true
		else -- check for the next part
			local a, b, c, d = parameters_cubic(graphs[1][1][1],
			graphs[1][1][2],graphs[1][math.floor(l/3)][1],
			graphs[1][math.floor(l/3)][2],
			graphs[1][math.floor(2*l/3)][1],
			graphs[1][math.floor(2*l/3)][2],
			graphs[1][l][1],graphs[1][l][2])
			return do_parameters_fit(a,b,c,d,"a*x^3+b*x^2+c*x+d",
			graphs[2],maxerror,false)
		end
	else
		return false
	end
end

-- returns true iff the inverse function is of type 
-- f(x)=a*x^3+b*x^2+c*x+d
-- a, b, c, d being real numbers
function is_cuberoot(graph,maxerror)
	local l = #graph
	local a, b, c, d = parameters_cubic(graph[1][2],graph[1][1],
	graph[math.floor(l/3)][2],graph[math.floor(l/3)][1],
	graph[math.floor(2*l/3)][2],graph[math.floor(2*l/3)][1],
	graph[l][2],graph[l][1])
	return do_parameters_fit(a,b,c,d,"a*x^3+b*x^2+c*x+d",graph,
	maxerror,true)
end

-- returns true iff the function is of type f(x)=a*x^3+b*x^2+c*x+d
-- a, b, c, d being real numbers
-- this takes several graph parts
-- the idea is to have a possibility to avoid tan(x)
function are_cuberoot(graphs,maxerror)
	if is_cuberoot(graphs[1],maxerror) then
		if #graphs < 2 then
			return true
		else -- check for the next part
			local a, b, c, d = parameters_cubic(graphs[1][1][2],
			graphs[1][1][1],graphs[1][math.floor(l/3)][2],
			graphs[1][math.floor(l/3)][1],
			graphs[1][math.floor(2*l/3)][2],
			graphs[1][math.floor(2*l/3)][1],
			graphs[1][l][2],graphs[1][l][1])
			return do_parameters_fit(a,b,c,d,"a*x^3+b*x^2+c*x+d",
			graphs[2],maxerror,true)
		end
	else
		return false
	end
end

-- returns true iff function is of type f(x)=a*x+b
-- a, b being real numbers
function is_affine(graph,maxerror)
	l = #graph
	local a, b = parameters_affine(graph[1][1],graph[1][2],
	graph[l][1],graph[l][2])
	return do_parameters_fit(a,b,0,0,"a*x+b",graph,maxerror,false)
end

-- what is the sum of the squared error
-- when comparing the bezier path
-- p.. control q and r .. s
-- with the graph g from index starti to endi
-- (looking at the points at roughly t=.33 and t=.67)
function squareerror(f,g,starti,endi,qx,qy,rx,ry)
	local result = 0
	for t = .33, .7, .34 do
		x = (1-t)^3*g[starti][1]+3*t*(1-t)^2*qx+3*t^2*(1-t)*rx+t^3*g[endi][1]
		y = (1-t)^3*g[starti][2]+3*t*(1-t)^2*qy+3*t^2*(1-t)*ry+t^3*g[endi][2]
		result = result + (y-f(x))^2
	end
	return result
end

function pointstobezier(qx,qy,rx,ry,sx,sy,rndx,rndy)
	return " .. controls (" .. round(qx,rndx) .. "," 
	.. round(qy,rndy) ..") and ("
	.. round(rx,rndx) .. "," 
	.. round(ry,rndy) .. ") .. (" 
	.. round(sx,rndx) .. "," 
	.. round(sy,rndy)..")"
end

-- take end points of a graph g of the function f
-- (from indices starti to endi)
-- without extrema or inflection points inbetween 
-- and try to approximate it with a cubic bezier curve
-- (round to rndx and rndy when printing)
function graphtobezierapprox(f,g,starti,endi,rndx,rndy,maxerror)
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
	local qx = px+.05*(cx-px)
	local qy = py+.05*(cy-py)
	local rx = sx+.05*(cx-sx)
	local ry = sy+.05*(cy-sy)
	local err = squareerror(f,g,starti,endi,qx,qy,rx,ry)
	for i = 2, 19 do
		for j = 2, 19 do
			xa = px+i*.05*(cx-px)
			ya = py+i*.05*(cy-py)
			xb = sx+j*.05*(cx-sx)
			yb = sy+j*.05*(cy-sy)
			-- now check, if xa and xb fit better
			-- at roughly t=0.33 and t=0.66 for f(x)
			-- than the last qx and rx did
			-- (sum of squares must be smaller)
			if squareerror(f,g,starti,endi,xa,ya,xb,yb) < err then
				qx = xa
				qy = ya
				rx = xb
				ry = yb
				err = squareerror(f,g,starti,endi,qx,qy,rx,ry)
			end
		end
	end
	-- check if it is close enough: (recycling err, xa, ya)
	err = 0
	for t = .1, .9, .1 do
		xa = (1-t)^3*g[starti][1]+3*t*(1-t)^2*qx+3*t^2*(1-t)*rx+t^3*g[endi][1]
		ya = (1-t)^3*g[starti][2]+3*t*(1-t)^2*qy+3*t^2*(1-t)*ry+t^3*g[endi][2]
		if abs(ya-f(xa)) > err then
			err = abs(ya-f(xa))
		end
	end
	if err <= maxerror then
		return pointstobezier(qx,qy,rx,ry,sx,sy,rndx,rndy)
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
		return graphtobezierapprox(f,g,starti,interindex,rndx,rndy,maxerror)
		.. graphtobezierapprox(f,g,interindex,endi,rndx,rndy,maxerror)
	end
end

-- like above but exact for quadratic and cubic (if not inverse)
-- resp. exact for squareroot and cuberoot (if inverse)
function graphtobezier(g,starti,endi,rndx,rndy,isinverse)
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
		return pointstobezier(qy,qx,ry,rx,sy,sx,rndy,rndx)
	else
		return pointstobezier(qx,qy,rx,ry,sx,sy,rndx,rndy)
	end
end

-- reverses a path p e.g. when p = "(0,1) -- (2,3)"
-- it returns "(2,3) -- (0,1)"
-- or when p = "(0,1) .. controls (2,3) and (4,5) .. (6,7)"
-- it returns "(6,7) .. controls (4,5) and (2,3) .. (0,1)"
function reversepath(p)
	local r = "" -- will become the reverse path
	local temppoint  ="" -- will store temporal points like "(0,1)"
	local tempbetween = "" -- will store things like " .. controls "
	for i = 1, #p do
		local c = string.sub(p,i,i)
		if c == "(" then
			if tempbetween == " .. " then
				r = " .. controls " .. r
			elseif tempbetween == " .. controls " then
				r = " .. " .. r
			else
				r = tempbetween .. r
			end
			tempbetween = ""
			temppoint = "("
		elseif c == ")" then
			r = temppoint .. ")" .. r
			temppoint = ""
		else
			if temppoint == "" then -- not reading a point
				tempbetween = tempbetween .. c
			else
				temppoint = temppoint .. c
			end
		end
	end
	return r
end

-- main function
function bezierplot(functionstring,xmin,xmax,ymin,ymax)
	local fstringreplaced = string.gsub(functionstring, "%*%*", "^")
	local f = assert(load("local x = ...; return " .. fstringreplaced)) 
	local isreverse = false
	if xmin > xmax then
		isreverse = true
	end
	xmin, xmax = math.min(xmin,xmax), math.max(xmin,xmax)
	local xstep = (xmax-xmin)/20000
	-- the output of the x coordinates will be rounded to rndx digits
	local rndx = math.max(0,math.floor(4.5-log(xmax-xmin)/log(10)))
	local xerror = abs(xmax-xmin)/(100*10^rndx)
	ymin, ymax = math.min(ymin,ymax), math.max(ymin,ymax)
	-- the output of the x coordinates will be rounded to rndy digits
	local rndy = math.max(0,math.floor(4.5-log(ymax-ymin)/log(10)))
	local yerror = (ymax-ymin)/(100*10^rndy)
	-- determine parts of the graph that are inside window
	local graphs = {}
	local outside = true -- value is outside window
	local i = 0
	local j = 0
	for n = 0, 20000 do
		local x = xmin + n/20000*(xmax-xmin)
		local y = f(x)
		if y >= ymin-yerror and y <= ymax+yerror then -- inside
			if outside then -- if it was outside before
				outside = false
				j = 0
				i = i + 1
				graphs[i] = {}
			end
			j = j + 1
			graphs[i][j] = {x,y}
		else
			outside = true
		end
	end

	local functiontype = "unknown"
	local bezierstring = ""

	-- go through the connected parts
	for part = 1, #graphs do 
		local d = diffgraph(f,graphs[part],xstep)
		-- check for type of function (only for the first part)
		if part == 1 then
			if is_affine(d,yerror) then
				functiontype = "affine"
			elseif are_cubic(graphs,yerror) then
				functiontype = "cubic"
			elseif are_cuberoot(graphs,xerror) then
				functiontype = "cuberoot"
			end
		end
		if functiontype ~= "cuberoot" then -- start with initial point
			bezierstring = bezierstring .. "(" .. round(d[1][1],rndx) 
			.. "," 	.. round(d[1][2],rndy) .. ")"
		end
		if functiontype == "affine" then 
			bezierstring = bezierstring .. " -- (" .. round(d[#d][1],
			rndx) .. "," .. round(d[#d][2],rndy) ..")"
		elseif functiontype == "cubic" then 
			local startindex = 1
			local extremainbetween = false
			for k = 2, #d do
				if d[k][5] then -- extrema 
					extremainbetween = true
					bezierstring = bezierstring 
					.. graphtobezier(d,startindex,k,rndx,rndy,false)
					startindex = k
				end
			end
			if not extremainbetween then
				for k = 2, #d do
					if d[k][6] then -- inflection point
						-- check, if the controlpoints are outside
						-- of the bounding box defined by the vertices
						-- (d[1][1],d[1][2]) and (d[#d][1],d[#d][2])
						local qx = d[1][1]+(d[#d][1]-d[1][1])/3
						local rx = d[1][1]+2*(d[#d][1]-d[1][1])/3
						local qy = d[1][2]+(qx-d[1][1])*d[1][3]
						local ry = d[#d][2]+(rx-d[#d][1])*d[#d][3]
						if math.max(qy,ry) > ymax 
						or math.min(qy,ry) < ymin then
							bezierstring = bezierstring ..graphtobezier(
							d,startindex,k,rndx,rndy,false)
							startindex = k
						end
					end
				end
			end
			if startindex ~= #d then -- if no special points inbetween
				bezierstring = bezierstring 
				.. graphtobezier(d,startindex,#d,rndx,rndy,false)
			end
		elseif functiontype == "cuberoot" then 
			-- we determine a, b, c, d and then
			-- get x' = 3ay^2+2by+c
			local a, b, c, dd = parameters_cubic(
			d[math.floor(.2*l)][2], d[math.floor(.2*l)][1],
			d[math.floor(.4*l)][2], d[math.floor(.4*l)][1],
			d[math.floor(.6*l)][2], d[math.floor(.6*l)][1],
			d[math.floor(.8*l)][2], d[math.floor(.8*l)][1])
			-- now recalculate the graph with the inverse function:
			-- we can increase the accuracy
			xstep = (ymax-ymin)/100000 -- inverse redefinition
			local finverse = assert(load("local x = ...; return "
			..a.."*x^3+"..b.."*x^2+"..c.."*x+"..dd))
			local graphinverse = {}
			local i = 1
			for y = ymin, ymax, xstep do
				local x = finverse(y)
				if x > xmin and x < xmax -- inside
				and abs(y-f(x)) < (ymax-ymin)/(100*10^rndy) then 
					graphinverse[i] = {y,x}
					i = i + 1
				end
			end
			d = diffgraph(finverse,graphinverse,xstep)
			bezierstring = bezierstring .. "(" .. round(d[1][2],rndy) 
			.. "," .. round(d[1][1],rndx) .. ")" -- initial point
			local startindex = 1
			for k = 2, #d do
				if d[k][6] then -- inflection point
					-- check, if the controlpoints are outside
					-- of the bounding box defined by the vertices
					-- (d[1][1],d[1][2]) and (d[#d][1],d[#d][2])
					local qx = d[1][1]+(d[#d][1]-d[1][1])/3
					local rx = d[1][1]+2*(d[#d][1]-d[1][1])/3
					local qy = d[1][2]+(qx-d[1][1])*d[1][3]
					local ry = d[#d][2]+(rx-d[#d][1])*d[#d][3]
					if math.max(qy,ry) > xmax 
					or math.min(qy,ry) < xmin then
						bezierstring = bezierstring..graphtobezier(
						d,startindex,k,rndx,rndy,true)
						startindex = k
					end
				end
			end
			if startindex ~= #d then -- if no special points inbetween
				bezierstring = bezierstring 
				.. graphtobezier(d,startindex,#d,rndx,rndy,true)
			end
		else	
			-- standard case (nothing special)			
			local startindex = 1
			for k = 2, #d do
				if d[k][5] or d[k][6] then -- extrema and inflection points
					bezierstring = bezierstring .. graphtobezierapprox(
					f,d,startindex,k,rndx,rndy,(ymax-ymin)/(0.5*10^rndy))
					startindex = k
				end
			end
			if startindex ~= #d then -- if no special points inbetween
				bezierstring = bezierstring .. graphtobezierapprox(f,d,
				startindex,#d,rndx,rndy,(ymax-ymin)/(0.5*10^rndy))
			end
		end
	end
	if isreverse then
		return reversepath(bezierstring)
	else
		return bezierstring
	end
end

-- main program --

if not pcall(debug.getlocal, 4, 1) then
	if #arg >= 1 then
		local xmin = -5
		local xmax = 5
		if #arg >= 2 then xmin = arg[2] end
		if #arg >= 3 then
			if arg[3] == arg[2] then
				xmax = xmin + 10
			else
				xmax = arg[3]
			end
		end
		local ymin = -5
		local ymax = 5
		if #arg >= 4 then ymin = arg[4] end
		if #arg >= 5 then 
			if arg[5] == arg[4] then
				ymax = ymin + 10
			else
				ymax = arg[5]
			end
		end
		print(bezierplot(arg[1],xmin,xmax,ymin,ymax))
	end
end



