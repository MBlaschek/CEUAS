import os
from datetime import datetime

SCP_SCRIPT_TEMPLATE = """#!/bin/bash
ssh {ssh_connect} 'mkdir {remote_dir}'
scp {local_dir}*.nc {ssh_connect}:{remote_dir}
scp {config} {ssh_connect}:{remote_home}/cdsobs_config.yml
scp {ingestion_script} {ssh_connect}:{remote_dir}
ssh {ssh_connect} 'chmod +x {remote_dir}ingestion_script.sh; {remote_dir}ingestion_script.sh'
"""

####       ####
# USER SETUP: #
#             #
ceuas_dir = '/srvfs/home/uvoggenberger/CEUAS/CEUAS/'
local_dir = "/mnt/users/scratch/uvoggenberger/CUON_HARVEST/202502/resort/2025/"
ssh_connect = "-J proxy@136.156.142.19  obs@dbdataset-cci2-0000"
remote_dir = "/mnt/public/" + local_dir.split("CUON_HARVEST")[-1].split("/resort")[0] + "/"
remote_home = "/home/obs"
#            #
##############  


def create_script(script_content, script_path):
    with open(script_path, 'w') as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)

def main():
    date = datetime.now().strftime("%Y%m%d")
    ingestion_script = ceuas_dir + "public/nrt_pipeline/ingestion_script.sh"


    config_file = f"{ceuas_dir}/public/nrt_pipeline/cdsobs_config_upload.yml"
    with open(f"{ceuas_dir}/public/nrt_pipeline/cdsobs_config.yml", 'r') as file:
        content = file.readlines()
    # Modify the content as needed
    modified_content = []
    for line in content:
        if "input_dir: " in line:
            line = f'input_dir: "{remote_dir}"\n'
        modified_content.append(line)
    # Write the modified content to a new file
    with open(config_file, 'w') as file:
        file.writelines(modified_content)

    # Paths for the generated scripts
    scp_script_path = "scp_upload.sh"

    # Create the SCP script
    scp_script = SCP_SCRIPT_TEMPLATE.format(
        ceuas_dir=ceuas_dir,
        ingestion_script=ingestion_script,
        ssh_connect=ssh_connect,
        local_dir=local_dir,
        remote_dir=remote_dir,
        remote_home=remote_home,
        config=config_file,
    )
    create_script(scp_script, scp_script_path)

    # Run the scripts
    # os.system(f"{scp_script_path}")

if __name__ == "__main__":
    main()