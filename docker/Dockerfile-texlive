## Start with the texlive image
FROM cboettig/texlive
MAINTAINER Matt Jones jones@nceas.ucsb.edu

# Create a volume containing texlive that other containers can mount
VOLUME /usr/local/texlive

# Add the texlive path from our linked container
ENV PATH $PATH:/usr/local/texlive/2014/bin/x86_64-linux/

# Configure supervisor to keep the container running
COPY docker/supervisord.conf /etc/supervisor/conf.d/supervisord.conf
RUN mkdir -p /var/log/supervisor \
  && chgrp staff /var/log/supervisor \
  && chmod g+w /var/log/supervisor \
  && chgrp staff /etc/supervisor/conf.d/supervisord.conf

CMD ["/usr/bin/supervisord", "-c", "/etc/supervisor/conf.d/supervisord.conf"]
