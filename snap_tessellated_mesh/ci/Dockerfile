ARG BASE=ubuntu:24.04
FROM $BASE

ARG N_PROCS=8

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get upgrade -y && apt-get install -y \
      gcc \
      gfortran \
      build-essential \
      wget \
      curl \
      libcurl4-gnutls-dev \
      bison \
      python3-dev \
      git \
      lcov \
      zlib1g-dev \
      libopenblas-dev \
      libdw-dev \
      libpfm4-dev \
      libunwind-dev \
      libgl1-mesa-dev \
      libgmp-dev \
      libxi-dev \
      libxmu-dev \
      tk-dev \
      tcl-dev \
      libfreetype6-dev \
      && \
      apt-get clean && rm -rf /var/lib/apt/lists/*

ENV PREFIX=/scratch
ENV ARCHIVE_DIR=${PREFIX}/archive
ENV SOURCE_DIR=${PREFIX}/source
ENV BUILD_DIR=${PREFIX}/build
ENV INSTALL_DIR=/opt

RUN mkdir -p ${PREFIX} && \
    cd ${PREFIX} && \
    mkdir archive && \
    mkdir source && \
    mkdir build

# Install CMake
RUN export CMAKE_VERSION=3.27.8 && \
    export  CMAKE_SHA256=dfedc30abe69dcabe326a4de632210bb52004185bd2c0b096c21924f4510681f && \
    export CMAKE_URL=https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-linux-x86_64.tar.gz && \
    export CMAKE_ARCHIVE=${ARCHIVE_DIR}/cmake.tar.gz && \
    export CMAKE_BUILD_DIR=${BUILD_DIR}/cmake && \
    wget --quiet ${CMAKE_URL} --output-document=${CMAKE_ARCHIVE} && \
    echo "${CMAKE_SHA256} ${CMAKE_ARCHIVE}" | sha256sum -c && \
    mkdir -p ${CMAKE_BUILD_DIR} && \
    tar xf ${CMAKE_ARCHIVE} -C ${CMAKE_BUILD_DIR} --strip-components=1 && \
    mv ${CMAKE_BUILD_DIR} ${INSTALL_DIR} && \
    rm -rf ${CMAKE_ARCHIVE} && \
    rm -rf ${CMAKE_BUILD_DIR}
ENV PATH=${INSTALL_DIR}/cmake/bin:$PATH

# Install OpenCascade 7.5.0
RUN export OPENCASCADE_HASH=29cf1a701ca269df10f2d1c4037e2826cf704355 && \
    export OPENCASCADE_URL=https://github.com/tpaviot/oce/archive/${OPENCASCADE_HASH}.tar.gz && \
    export OPENCASCADE_ARCHIVE=${ARCHIVE_DIR}/opencascade.tar.gz && \
    export OPENCASCADE_SOURCE_DIR=${SOURCE_DIR}/opencascade && \
    export OPENCASCADE_BUILD_DIR=${BUILD_DIR}/opencascade && \
    export OPENCASCADE_INSTALL_DIR=${INSTALL_DIR}/opencascade && \
    wget --quiet ${OPENCASCADE_URL} --output-document=${OPENCASCADE_ARCHIVE} && \
    mkdir -p ${OPENCASCADE_SOURCE_DIR} && \
    tar -xf ${OPENCASCADE_ARCHIVE} -C ${OPENCASCADE_SOURCE_DIR} --strip-components=1 && \
    mkdir -p ${OPENCASCADE_INSTALL_DIR} && \
    mkdir ${OPENCASCADE_BUILD_DIR} && \
    cd ${OPENCASCADE_BUILD_DIR} && \
    cmake -D OCE_INSTALL_PREFIX=${OPENCASCADE_INSTALL_DIR} \
          -D OCE_TESTING=OFF \
          -D OCE_VISUALISATION=OFF \
          -D OCE_DISABLE_X11=ON \
          ${OPENCASCADE_SOURCE_DIR} && \
    make -j${N_PROCS} install && \
    rm -rf ${OPENCASCADE_ARCHIVE} && \
    rm -rf ${OPENCASCADE_BUILD_DIR} && \
    rm -rf ${OPENCASCADE_SOURCE_DIR}
ENV OPENCASCADE_DIR=${INSTALL_DIR}/opencascade

# Install deal.II 9.5.1
RUN export DEAL_II_HASH=9e847302b21355f355c87890477f4dd485da26b1 && \
    export DEAL_II_URL=https://github.com/dealii/dealii/archive/${DEAL_II_HASH}.tar.gz && \
    export DEAL_II_ARCHIVE=${ARCHIVE_DIR}/dealii.tar.gz && \
    export DEAL_II_SOURCE_DIR=${SOURCE_DIR}/dealii && \
    export DEAL_II_BUILD_DIR=${BUILD_DIR}/dealii && \
    export DEAL_II_INSTALL_DIR=${INSTALL_DIR}/dealii && \
    wget --quiet ${DEAL_II_URL} --output-document=${DEAL_II_ARCHIVE} && \
    mkdir -p ${DEAL_II_SOURCE_DIR} && \
    tar -xf ${DEAL_II_ARCHIVE} -C ${DEAL_II_SOURCE_DIR} --strip-components=1 && \
    mkdir -p ${DEAL_II_BUILD_DIR} && cd ${DEAL_II_BUILD_DIR} && \
    cmake -D CMAKE_BUILD_TYPE=DebugRelease \
          -D CMAKE_CXX_STANDARD=17 \
          -D DEAL_II_ALLOW_AUTODETECTION=OFF \
          -D DEAL_II_WITH_COMPLEX_VALUES=OFF \
          -D DEAL_II_WITH_OPENCASCADE=ON \
          -D OPENCASCADE_DIR=${OPENCASCADE_DIR} \
          -D DEAL_II_COMPONENT_EXAMPLES=OFF \
          -D CMAKE_INSTALL_PREFIX=${DEAL_II_INSTALL_DIR} \
          ${DEAL_II_SOURCE_DIR}  && \
    make -j${N_PROCS} install && \
    rm -rf ${DEAL_II_ARCHIVE} && \
    rm -rf ${DEAL_II_BUILD_DIR}
# We keep the source file for debugging purpose
ENV DEAL_II_DIR=${INSTALL_DIR}/dealii
