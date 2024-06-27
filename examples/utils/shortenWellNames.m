function schedule = shortenWellNames(schedule, f)
    % Utility to shorten well names
    % f is a function handle s.t. newname = f(oldname)
    for i = 1:numel(schedule.control(1).W)
        for j = 1:numel(schedule.control)
            schedule.control(j).W(i).name = f(schedule.control(j).W(i).name);
        end
    end 
end
