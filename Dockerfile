## Start with the ropensci image, and put our check code inline
FROM rocker/ropensci
MAINTAINER Matt Jones jones@nceas.ucsb.edu

# Install packages needed for R checks
RUN apt-get update \
  && apt-get install -y --no-install-recommends -t unstable \
    qpdf

# Copy check script into the container
COPY docker /home/rstudio
RUN chown -R rstudio /home/rstudio

#ENTRYPOINT ["R", "CMD", "check", "--as-cran", "codyn"]
#CMD bash -c /home/rstudio/codyn-check.sh
CMD cd /home/rstudio && chmod ugo+rx c* && cp Rprofile /root/.Rprofile && ls -l && id && cd /src/dev && ls -l && R CMD build codyn && R CMD check --as-cran codyn_0.9.0.9000.tar.gz
#CMD ['/bin/sh']

RUN echo "NEXT STEPS:" && \
    echo "Build & run the container to execute the tests:" && \
    echo "$ docker build -t metamattj/codyn-check ." && \
    echo '$ docker run --rm -v $(pwd):/src/dev/codyn --name=codyn -it metamattj/codyn-check'
