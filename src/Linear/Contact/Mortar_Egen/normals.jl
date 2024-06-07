function calculate_normals(elements, X)
    normals = empty(X)
    for (elid, elcon) in elements
        d = X[elcon[2]] - X[elcon[1]]
        # Riktning beror på hur elementnoder är orienterade
        if d[1] < 0
        	d = -d
        end
        n = [-d[2], d[1]]
        n /= norm(n) # Här borde vi normalen reversas om nödvändigt, kan använda elid = "element id"
        for nid in elcon
            haskey(normals, nid) || (normals[nid] = zeros(2))
            normals[nid] += n
        end
    end
    for (nid, n) in normals
        normals[nid] /= norm(normals[nid])
    end
    return normals
end
