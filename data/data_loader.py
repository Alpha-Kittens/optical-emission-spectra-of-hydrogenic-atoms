import numpy as np

def read_data(file_path):
    """
    Reads scanning data from given file.

    Arguments: 
        * `file_path` (str): string containing the relative location of the file from the directory containing the function call.

    Returns: 
        * `data` (nx2 numpy array): numpy array (matrix) with 2 columns:
            - First column is wavelength in Angstroms
            - Second column is measured counts per second
    """
    
    data = np.array([[0,0]])

    #print ("starting to read")
    
    with open(file_path) as f:
        
        for row in f.readlines():

            line = row.split("\t")

            #print (line)

            data = np.append(data, [[float(line[1]), float(line[2])]], axis = 0)

    return data[1:]

if __name__ == "__main__":

    filepath = "./day_1/calibration_scan_1_3800_4600"

    print(read_data(filepath))

