function penalty(g::Float64)
    return ϵ * max(0,g)
end


# Skriv på ett sätt så att slave och master nodes bara behöver tillhöra olika set, 
# då kan självkontakt inkluderas vid ett senare skede....
function contact_search(slave_nodes,master_nodes,slave_elems,master_elems)
    Pairs = Dict;


    return Pairs
end