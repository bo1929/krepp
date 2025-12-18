FROM alpine:3.13.5
RUN apk update && \
    apk add bash g++ git make zlib-dev && \
    git clone --recurse-submodules https://github.com/bo1929/krepp.git && \
    cd krepp && \
    make && \
    mv krepp /usr/local/bin/krepp && \
    cd .. && \
    rm -rf krepp
