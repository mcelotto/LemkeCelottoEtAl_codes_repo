function index = median_time(time_profile)

index = find(cumsum(time_profile) >= sum(time_profile) / 2, 1, 'first');

end