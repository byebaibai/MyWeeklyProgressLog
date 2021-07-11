resultMatrix = ScriptA(32, 3);
file = fopen("ScriptC.txt", "w");
colNames = ["col1" "col2" "col3"];
for ci = 1:size(resultMatrix, 2)
    fprintf(file, "%s\t", colNames(ci));
end
fprintf(file, "\n");
for ri = 1:size(resultMatrix, 1)
    for ci = 1:size(resultMatrix, 2)
        fprintf(file, "%d\t", resultMatrix(ri, ci));
    end
    fprintf(file, "\n");
end