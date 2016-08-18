FROM ubuntu:14.04

MAINTAINER eric.t.dawson@gmail.com

USER root
RUN apt-get -yq update && \
    apt-get -yq install build-essential libzlib gcc4.9

RUN make
