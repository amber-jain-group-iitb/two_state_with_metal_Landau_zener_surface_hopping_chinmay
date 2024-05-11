
import fileinput

inputfil="AFSSH.inp"
with fileinput.FileInput(inputfil, inplace=True) as file:
    for line in file:
        print(line.replace("2 !! iflow", "3 !! iflow"), end='')
