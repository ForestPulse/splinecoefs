FROM ghcr.io/osgeo/gdal:ubuntu-small-3.12.2 AS builder

# disable interactive frontends
ENV DEBIAN_FRONTEND=noninteractive 

# Environment variables
ENV SOURCE_DIR=/source_dir

# Copy src to SOURCE_DIR
RUN mkdir -p $SOURCE_DIR/src
WORKDIR $SOURCE_DIR
COPY . $SOURCE_DIR/src

RUN apt-get -y update && apt-get -y upgrade && \
  apt-get -y install \
    build-essential \
    wget && \
  apt-get clean && rm -r /var/cache/

# build GSL from source
RUN mkdir -p $SOURCE_DIR/gsl && cd $SOURCE_DIR/gsl && \
wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.8.tar.gz && \
tar -xvf gsl-2.8.tar.gz && \
cd gsl-2.8 && \
./configure --prefix=/opt/libgsl28 && \
make -j7 && \
make install && \
ldconfig


# Build, install
RUN echo "building splinecoefs" && \
  cd $SOURCE_DIR/src && \
  make



FROM ghcr.io/osgeo/gdal:ubuntu-small-3.12.2 AS final

# Environment variables
ENV SOURCE_DIR=/source_dir

COPY --from=builder /opt/libgsl28 /opt/libgsl28
COPY --from=builder $SOURCE_DIR/src/splinecoefs /usr/local/bin


# Create a dedicated 'docker' group and user
RUN groupadd docker && \
useradd -m docker -g docker -p docker && \
chmod 0777 /home/docker && \
chgrp docker /usr/local/bin && \
mkdir -p /home/docker/bin && chown docker /home/docker/bin

# Use this user by default
USER docker
ENV HOME=/home/docker
ENV PATH="$PATH:/home/docker/bin"
WORKDIR /home/docker

CMD ["splinecoefs"]
