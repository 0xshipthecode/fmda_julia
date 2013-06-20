using Stations
s = Stations.load_station_info("espc2.info")
Stations.load_station_data(s, "espc2_1.csv")

t = Stations.obs_times(s)
o = Stations.observations(s, "RH")
o1 = o[t[1]]

println(o1)
