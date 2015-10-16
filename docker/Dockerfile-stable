## Start with the ropensci image, and put our check code inline
FROM rocker/ropensci
MAINTAINER Matt Jones jones@nceas.ucsb.edu

# Set up the working directory to a locally mounted volume
WORKDIR /src

# Install packages needed for R checks
RUN apt-get update \
  && apt-get install -y --no-install-recommends -t unstable \
    qpdf

# Copy check script into the container
COPY docker /root

CMD pwd
