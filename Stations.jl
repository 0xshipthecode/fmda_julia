module Stations

using Utils
import Utils.GeoLoc

using Calendar
import Calendar.CalendarTime

import Base.show


abstract AbstractStation


type Observation

    # station of origin
    station::AbstractStation

    # observation time
    tm::CalendarTime

    # observed value
    value::Float64

    # name of the variable
    obs_type::String

    # variance of the observation
    var :: Float64

    # default constructor
    Observation(s, tm, value, obs_type, var) = new(s, tm, value, obs_type, var)

end

type Station <: AbstractStation

    # station id
    id::String

    # station name
    name::String

    # station location (lat/lon decimal)
    loc::GeoLoc

    # nearest grid point (tuple)
    ngp::(Int64, Int64)

    # station elevation (meters)
    elevation::Float64

    # dictionary mapping variable names to observations
    obs::Dict{String,Array{Observation}}

    # available observation types from the station
    sensors :: Array{String}

    Station() = new("", "", GeoLoc(0.0, 0.0), (-1,-1), 0.0,  Dict{String,Array{Observation}}(), Array(String,0))

    Station(id, name, loc, elevation) = new(id, name, loc, (-1,-1), elevation, Dict{String,Array{Observation}}(), Array(String,0))

end

obs_variance(o::Observation) = o.var
obs_value(o::Observation) = o.value
obs_type(o::Observation) = o.obs_type
obs_station(o::Observation) = o.station
obs_station_id(o::Observation) = o.station.id
nearest_grid_point(o::Observation) = o.station.ngp


function show(io::IO, o::Observation)
#    write(io, "Obs($o.value of $o.obs_type at $o.station.id, $o.tm with var $o.var)")
    print(io, "Obs(", o.obs_type, "=", o.value, " at ", o.station.id, ", ", o.tm, " with var ", o.var, ")")
end


obs_times(s::Station, obs_type::String) = [o.tm for o in observations(s, obs_type)]
obs_var(s::Station, obs_type::String) = s.obs[obs_type]
id(s::Station) = s.id
name(s::Station) = s.name
observations(s::Station, obs_type::String) = s.obs[obs_type]

function register_to_grid(s::Station, lat::Array{Float64,2}, lon::Array{Float64,2})
    s.ngp = ind2sub(size(lat), indmin( (s.loc.lat - lat).^2 + (s.loc.lon - lon).^2 ))
end
    


function observations(s::Station, obs_type::String, tm::Array{CalendarTime})
    
    # retrieve all the observations for the given type
    obs = s.obs[obs_type]
    otm = [o.tm for o in obs]
    obs_ret = Dict(CalendarTime, Observation)
    
    i, j = 1, 1
    while (i <= length(tm)) && (j <= length(otm))
        if tm[i] == otm[j]
            obs_ret[tm[i]] = obs[i]
            i += 1
            j += 1
        elseif tm[i] > obst[j]
            j += 1
        elseif tm[i] < obst[j]
            j += 1
        end
    end
    
    return obs_ret
end



function load_station_info(fname :: String)

    s = Station()

    # open file
    io = open(fname, "r")

    # read station id
    s.id = strip(readline_skip_comments(io))

    # read station name
    s.name = strip(readline_skip_comments(io))

    # read station geo location
    s.loc = Utils.parse_geolocation(readline_skip_comments(io))

    # read elevation
    s.elevation = float(strip(readline_skip_comments(io)))

    # read all observation types acquired by the station
    s.sensors = map(x -> strip(x), split(readline_skip_comments(io), ","))
    println(s.sensors)

    # create var time series containers
    for v in s.sensors
        s.obs[v] = Array(Float64, 0)
    end

    # close the info file
    close(io)

    return s
end



function load_station_data(s::Station, fname::String)
    # open file
    io = open(fname, "r")

    # read observations, time is first (and in GMT), then one value per measurement or nan if not available
    while !eof(io)

        # read the date (if cannot read date, file is empty)
        tm_line = readline_skip_comments(io)
        if tm_line == nothing  break end
        tm_str = strip(tm_line)
        tm = Calendar.parse("yyyy-MM-dd_HH:mm", tm_str, "GMT")

        # read in observed variables
        variables = map(x -> strip(x), split(readline_skip_comments(io), ","))

        # read in observed values
        vals = map(x -> float(x), split(readline_skip_comments(io), ","))

        # read in variances of observations
        variances = map(x -> float(x), split(readline_skip_comments(io), ","))

        # add each observation to the station dictionary
        for i in 1:length(variables)
            var = variables[i]
            push!(s.obs[var], Observation(s, tm, vals[i], var, variances[i]))
        end

    end

    # cleanup
    close(io)

    # return the created station
    return s

end


function build_observation_data(ss::Array{Station}, obs_type::String)
    Ns = length(ss)
    obs = Observation[]
    
    # observation data
    for s in ss
        append!(obs, observations(s, obs_type))
    end

    # repackage into a time-indexed structure
    obs_data = Dict{CalendarTime,Array{Observation}}()
    for o in obs
        if has(obs_data, o.tm) push!(obs_data[o.tm], o) else obs_data[o.tm] = [o] end
    end

    return obs_data
end


function readline_skip_comments(io :: IO)
    s = "#"
    while (length(s) == 0) || (s[1] == '#')
        if !eof(io)
            s = readline(io)
        else
            return nothing
        end
    end
    return s
end


end
