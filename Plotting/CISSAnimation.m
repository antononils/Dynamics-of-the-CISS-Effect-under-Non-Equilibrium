function CISSAnimation(t, sites, spinUpData, spinDownData, title)
    % Function for animating probability densities over time
    
    % Create a figure
    figure

    % Loop through time and plot probability density for each spin
    for k = 1:numel(t)-1
        % Plot the data for current time
        plot(sites, spinDownData(k,:), sites, spinUpData(k, :),'-')
        
        % Create labels for data and axis
        legend("Spin down", "Spin up")
        xlabel("Site Index")
        ylabel("Probability Density")

        % Save frame for current time
        F(k) = getframe(gcf);

        % Pause time to achieve "video-like" animation
        pause(0.0001)
    end
    
    % Create object to save video in
    writerObj = VideoWriter(title);
    writerObj.FrameRate = 60;
    open(writerObj);
    
    % Loop through frames and convert to an animation
    for i=1:length(F)
        frame = F(i);    
        writeVideo(writerObj, frame);
    end
    close(writerObj)
end