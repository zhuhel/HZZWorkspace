# to build an image locally,
# go into HZZWorkspace (top level)
# and call docker build .

FROM atlas/analysisbase:latest

ADD . /hzz/source/HZZWorkspace
WORKDIR /hzz/build

RUN source ~/release_setup.sh &&  \
    sudo chown -R atlas /hzz && \
    mv ../source/HZZWorkspace/CMakeLists.topLevel.txt ../source/CMakeLists.txt && \
    cmake ../source && \
    make -j4 && \
    echo "source /hzz/source/HZZWorkspace/SetupInDocker.sh" >> /home/atlas/.bashrc
