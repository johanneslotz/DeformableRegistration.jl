function *(factor,JTuple)
    return factor.*JTuple[1],factor.*JTuple[2],factor.*JTuple[3]
end

function +(JTuple1,JTuple2)
    return JTuple1[1].+JTuple2[1],JTuple1[2].+JTuple2[2],JTuple1[3].+JTuple2[3]
end
