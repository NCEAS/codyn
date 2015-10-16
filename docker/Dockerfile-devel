## Start with the rstudio daily image, which is built on a daily R build snapshot
FROM rocker/rstudio-daily:verse
MAINTAINER Matt Jones jones@nceas.ucsb.edu

# Add the texlive path from our linked container
ENV PATH $PATH:/usr/local/texlive/2014/bin/x86_64-linux/

# Copy check script into the container
COPY docker /root

CMD pwd
