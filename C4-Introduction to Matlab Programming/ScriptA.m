%{
    Problem Intro.
1. Create a 4 x 8 matrix of randomly generated numbers.

2. Loop through all rows and columns, and test whether each element is greater than 0.5.

3. Report the results of the test along with the value of the matrix element and its 
row-column position. For example, your Matlab script should print The 3rd row and 8th 
column has a value of 0 .42345 and is not bigger than 0.5.

4. Make sure to add exceptions to print out 1st, 2nd, and 3rd, instead of lth, 2th, and 3th.

5. Put this code into a separate function that you can call from the command line with two inputs, 
corresponding to the number of rows and the number of columns of the matrix.
%}

function labelMatrix = ScriptA(rows, cols)
    randomMatrix = rand(rows, cols);
    labelMatrix = zeros(rows, cols);
    for rowi = 1:size(randomMatrix, 1)
        for coli = 1:size(randomMatrix, 2)
            if randomMatrix(rowi, coli) > 0.5
                fprintf("The %d%s row and %d%s column has a value of %.5f and is bigger than 0.5\n", ...
                rowi, getRankSuffix(rowi), coli, getRankSuffix(coli), randomMatrix(rowi, coli))

                labelMatrix(rowi, coli) = 1;
            end
        end
    end

    %getRankSuffix - Description
    %
    % Syntax: rankSuff = getRankSuffix(rank)
    %
    % get string suffix according to rank number
    function rankSuf = getRankSuffix(theRank)
        switch theRank
            case 1
                rankSuf = "st";
            case 2
                rankSuf = "nd";
            case 3
                rankSuf = "rd";
            otherwise
                rankSuf = "th";
        end
    end
end





