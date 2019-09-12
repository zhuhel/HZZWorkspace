# Tell people what they need to set up and point to documentation

echo "============================================="
echo "You can safely ignore the above; "
echo "this container is already fully set up"
echo "to start using the Workspace framework. "
echo
echo "For further documentation, please refer to  "
echo "      https://gitlab.cern.ch/HZZ/HZZSoftware/HZZWorkspace/wikis/home"
echo
echo "You are in the HZZ build directory."
echo "The code is mounted in the Docker volume as /hzz/"
echo "Thank you for using HZZWorkspace with Docker! "
echo "============================================="
source /home/atlas/release_setup.sh
source ${AnalysisBase_PLATFORM}/setup.sh
python ../source/HZZWorkspace/verifysetup.py
