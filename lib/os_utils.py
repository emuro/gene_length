import os  # to deal with the dirs and files



def recursiveSearch_leave_subdirs(root_path="./", path_list=[]):
    ''' given the next inputs:
            -root_path
        Find all the subdirectories from the path that are leaves. That is
        nothing in way up to the leave-dirs.

        output:
            add the path to path_list

        Note on computational time: takes some time ~2sec for around 3*10^final dirs (more dirs in the middle)
    '''

    n1_withpath=list()
    for dI in os.listdir(root_path):
        if os.path.isdir(os.path.join(root_path,dI)):
            n1_withpath.append(os.path.join(root_path,dI))
    if len(n1_withpath)==0:
        path_list.append(root_path)
        return
    else:
        for i in range(0, len(n1_withpath)):
            recursiveSearch_leave_subdirs(n1_withpath[i], path_list)