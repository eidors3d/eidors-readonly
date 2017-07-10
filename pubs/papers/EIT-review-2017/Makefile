all:
	make -C docs
	cp docs/EIT-review.pdf .
	make -C docs EIT-review.zip
	cp docs/EIT-review.zip .

clean:
	make -C docs clean
	-rm EIT-review.pdf
	-rm EIT-review.zip
