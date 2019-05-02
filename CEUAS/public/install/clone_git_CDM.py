""" Module for cloning the CDM GitHub repository                                                       
                                                                                                       
    Author      :: Ambrogi Federico , federico.ambrogi@univie.ac.at                                                                                                                                       
"""

import os,sys
import git


destination = "../common/CDM"

if not os.path.isdir(destination):
    os.mkdir(destination)

print("Cloning the CDM GitHub repository in %s" %destination)
git.Git("").clone("https://github.com/glamod/common_data_model.git")
print("Done!")
