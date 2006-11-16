# Build script for EIDORS Packages

# Need to have a local apache pointing to eidors

release:
	mkdir release


release/tutorial: release
	mkdir -p tmp && cd tmp && \
	wget -r http://localhost:8080/eidors/tutorial
	find tmp/localhost+8080/eidors/tutorial -name \*.shtml -exec \
	 perl -i -pe's{\.shtml\s*"}{.html"}' \{} \;
	find tmp/localhost+8080/eidors/tutorial -name \*.shtml -print \
         | while read a ; do echo $$a `basename $$a` ; done
