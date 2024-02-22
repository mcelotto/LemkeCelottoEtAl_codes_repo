function com = center_of_mass_time(time_profile)

t_points = 1:numel(time_profile);
com = sum(time_profile.*t_points)/sum(time_profile);

end