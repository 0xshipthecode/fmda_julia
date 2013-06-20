module Utils


type GeoLoc
    lat :: Float64
    lon :: Float64

    GeoLoc(lat,lon) = new(lat, lon)
end


function show(io :: IO, gl :: GeoLoc)
    print("GL($gl.lat, $gl.lon)")
end


function parse_geolocation(s::String)
    lat, lon = map(x -> float(strip(x)), split(s, ","))
    return GeoLoc(lat, lon)
end


end