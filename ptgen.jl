function halton(b, i)
        f = 1
        r = 0
        while i > 0
                f = f / b
                r = r + f* (mod(i, b))
                i = floor(i/b)
        end
        return r
end

function haltons(b1, b2, N)
        points = zeros(N,2)
        for i=1:N
                points[i, 1] = halton(b1, i)
                points[i, 2] = halton(b2, i)
        end
        return points
end

# Enters N halton point into domain xf with bounding box with right end vector v
function EnterHalton(b1, b2, xf, v, N)
        pts = zeros(N, 2)
        i = 1
        j = 1
        while i < N+1
                pt = v.*[halton(b1, j), halton(b2, j)]
                if xf(pt) > 0.5
                        pts[i, :] = pt
                        i = i+1
                end
                j = j + 1
        end
        return pts
end


# Enters NxK points on a regular grid into the 
# interior of a rectangle with bounding box right end
# vector v
function EnterEqui(v, N, K)
        pts = zeros(N*K, 2)
        
        for n=1:N
                for k=1:K
                        pts[(n-1)*K + k, 1] = n/(N+1)*v[1]
                        pts[(n-1)*K + k, 2] = k/(K+1)*v[2]
                end
        end

        return pts
end

