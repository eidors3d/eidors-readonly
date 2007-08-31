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

# TODO FOR RELEASE
find -name \*.shtml | while read a  ; do mv $a `dirname $a`/`basename $a .shtml`.html ; done
find -name \*.html -exec perl -i -pe's{\.shtml\s*"}{.html"}' \{} \;
find -name \*.html.bak | xargs rm
# FIND LINKS THAT ARE TOO BIG (>0.5MB)
#  cg_deforming_tank_phantom
#  cg_normal_breathing
#  netgen_moving_ball
