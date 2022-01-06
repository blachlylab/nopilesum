FROM charesgregory/dlang-htslib-static

ADD . /home/nopilesum

WORKDIR /home/nopilesum

RUN dub build --compiler ldc2 -c static-alpine -b release

RUN cp nopilesum /usr/local/bin

ENTRYPOINT ["/usr/local/bin/nopilesum"]