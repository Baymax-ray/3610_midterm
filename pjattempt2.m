function monteCarloSimulation(population, vaccinesPerMonth, deathRateByAge, infectionRateByAge)
    population = 107000; % Total population
    vaccinesPerMonth = 3500; % Number of vaccines to distribute per month
    deathRateByAge = [0.05, 0.02, 0.01]; % Death rate by age group
    infectionRateByAge = [0.1, 0.05, 0.02]; % Infection rate by age group

    numAgeGroups = numel(deathRateByAge);

    numMonths = ceil(population / vaccinesPerMonth);

    totalPeopleAffected = zeros(1, numMonths);
    totalPeopleKilled = zeros(1, numMonths);
    totalYearsOfLifeLost = zeros(1, numMonths);
    totalPeopleVaccinatedByAge = zeros(numAgeGroups, numMonths);

    cumulativePeopleAffected = 0;
    cumulativePeopleKilled = 0;
    cumulativeYearsOfLifeLost = 0;
    cumulativePeopleVaccinatedByAge = zeros(numAgeGroups, 1);

    for month = 1:numMonths
        % Calculate the number of vaccines available this month
        vaccinesAvailable = min(vaccinesPerMonth, population);

        % Determine the order of vaccination by age group (oldest to youngest)
        orderOfVaccination = 1:numAgeGroups;

        % Shuffle the order randomly
        orderOfVaccination = orderOfVaccination(randperm(numAgeGroups));

        % Initialize variables to track the impact
        peopleAffectedByAge = zeros(1, numAgeGroups);
        peopleKilledByAge = zeros(1, numAgeGroups);
        yearsOfLifeLostByAge = zeros(1, numAgeGroups);
        peopleVaccinatedByAge = zeros(1, numAgeGroups);

        % Simulate the disease spread and vaccination for each age group
        for ageGroup = orderOfVaccination
            % Calculate the number of susceptible people in this age group
            susceptible = floor(population / numAgeGroups);

            % Calculate the number of infected people based on infection rate
            infected = floor(susceptible * infectionRateByAge(ageGroup));

            % Calculate the number of people affected (can't be higher than susceptible)
            affected = min(infected, susceptible);

            % Calculate the number of people killed based on death rate
            killed = floor(affected * deathRateByAge(ageGroup));

            % Calculate years of life lost
            yearsOfLifeLost = killed * (80 - ageGroup); % Assuming a life expectancy of 80 years

            % Calculate the number of people to vaccinate, excluding those who are ill
            peopleToVaccinate = min(vaccinesAvailable, susceptible - affected);

            % Update totals for this month
            peopleAffectedByAge(ageGroup) = affected;
            peopleKilledByAge(ageGroup) = killed;
            yearsOfLifeLostByAge(ageGroup) = yearsOfLifeLost;
            peopleVaccinatedByAge(ageGroup) = peopleToVaccinate;

            % Update total population, vaccines available, and health status
            population = population - affected;
            vaccinesAvailable = vaccinesAvailable - peopleToVaccinate;
            % Remove ill people from population (they can't be vaccinated)
            population = population - affected;

            % Update cumulative totals
            cumulativePeopleAffected = cumulativePeopleAffected + affected;
            cumulativePeopleKilled = cumulativePeopleKilled + killed;
            cumulativeYearsOfLifeLost = cumulativeYearsOfLifeLost + yearsOfLifeLost;
            cumulativePeopleVaccinatedByAge(ageGroup) = cumulativePeopleVaccinatedByAge(ageGroup) + peopleToVaccinate;
        end

        % Store the cumulative results for this month
        totalPeopleAffected(month) = cumulativePeopleAffected;
        totalPeopleKilled(month) = cumulativePeopleKilled;
        totalYearsOfLifeLost(month) = cumulativeYearsOfLifeLost;
        totalPeopleVaccinatedByAge(:, month) = cumulativePeopleVaccinatedByAge;
    end

    % Plot cumulative results
    figure;

    subplot(4, 1, 1);
    plot(1:numMonths, totalPeopleAffected, 'b', 'LineWidth', 2);
    title('Cumulative Total Number of People Affected');
    xlabel('Months');
    ylabel('Count');
    
    subplot(4, 1, 2);
    plot(1:numMonths, totalPeopleKilled, 'r', 'LineWidth', 2);
    title('Cumulative Total Number of People Killed');
    xlabel('Months');
    ylabel('Count');
    
    subplot(4, 1, 3);
    plot(1:numMonths, totalYearsOfLifeLost, 'g', 'LineWidth', 2);
    title('Cumulative Total Years of Life Lost');
    xlabel('Months');
    ylabel('Years');
    
    subplot(4, 1, 4);
    plot(1:numMonths, totalPeopleVaccinatedByAge', 'LineWidth', 2);
    title('Cumulative Total People Vaccinated by Age Group');
    xlabel('Months');
    ylabel('Count');
    legend('Age Group 1', 'Age Group 2', 'Age Group 3');
end
