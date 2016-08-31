def write_tecplot_scatter(filename,title,variable_names,data_points):

    # Open the data file
    fid = open(filename,'w')
    
    # Write the title
    fid.write('title = '+title+'\n')

    # Write tha variable names
    varnames_commas = ','.join(variable_names) # Merge names in a single string separated by commas
    fid.write('variables = '+varnames_commas+',\n') # Write to the file

    # Write data points
    if type(data_points) is list: # Check if user provided a list
        for point in data_points:
            str_points = [str(x) for x in point] # Covert each entry to string
            str_points = ' '.join(str_points) # Merge all entries in a single string separated by whitespace
            fid.write(str_points+'\n') # Write to file
    else: # The user probably provided a numpy array
        for index in range(data_points.shape[0]):
            str_points = [str(x) for x in data_points[index,:]]
            str_points = ' '.join(str_points) # Merge all entries in a single string separated by whitespace
            fid.write(str_points+'\n') # Write to file
    
    # Close file
    fid.close()
