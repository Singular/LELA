/* linbox/util/splicer.C
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 */

#include <iostream>
#include <sstream>
#include <algorithm>

#include "linbox/util/splicer.h"
#include "linbox/util/commentator.h"

namespace LinBox
{

std::ostream &operator << (std::ostream &os, const Block &b)
{
	os << "Block (source=" << b.source () << ",source-index=" << b.sourceIndex () << ",dest-index=" << b.destIndex () << ",size=" << b.size () << ")";

	return os;
}

struct CompareBlockSources { inline bool operator () (const Block &b1, const Block &b2) { return b1.source () < b2.source (); } };

std::vector<Block> &Splicer::map_blocks (std::vector<Block> &output, const std::vector<Block> &outer_blocks, const std::vector<Block> &inner_blocks, unsigned int inner_source) const
{
	linbox_check (!inner_blocks->empty ());

	std::vector<Block>::const_iterator outer_block;
	std::vector<Block>::const_iterator inner_block = inner_blocks.begin ();

	size_t curr_dest_idx = 0, rest_size = 0, curr_size = 0;
	std::map<size_t, size_t> curr_source_idx;

	commentator.start ("Mapping blocks", __FUNCTION__);

	size_t max_inner_block = (std::max_element (inner_blocks.begin (), inner_blocks.end (), CompareBlockSources ()))->source ();

	for (outer_block = outer_blocks.begin (); outer_block != outer_blocks.end (); ++outer_block) {
		if (outer_block->source () != inner_source) {
			if (outer_block->source () < inner_source)
				output.push_back (*outer_block);
			else
				output.push_back (Block (outer_block->source () - inner_source + max_inner_block, outer_block->sourceIndex (), outer_block->destIndex (), outer_block->size ()));

			curr_dest_idx += outer_block->size ();

			commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION)
				<< "New block created (from outer block): " << output.back () << std::endl;
		} else {
			if (rest_size > 0) {
				curr_size = std::min (rest_size, outer_block->size ());

				output.push_back (Block (inner_block->source () + outer_block->source (), curr_source_idx[inner_block->source ()], curr_dest_idx, curr_size));

				commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION)
					<< "New block created (leftover): " << output.back () << std::endl
					<< "Inner source-block: " << *inner_block << std::endl
					<< "Outer source-block: " << *outer_block << std::endl;

				curr_source_idx[inner_block->source ()] += curr_size;
				curr_dest_idx += curr_size;
				rest_size -= curr_size;

				if (rest_size == 0)
					++inner_block;
			}

			while (inner_block != inner_blocks.end () && rest_size == 0 && outer_block->destIndex () + outer_block->size () > curr_dest_idx) {
				curr_size = std::min (inner_block->size (), outer_block->destIndex () + outer_block->size () - curr_dest_idx);

				if (curr_source_idx.find (inner_block->source ()) == curr_source_idx.end ())
					curr_source_idx[inner_block->source ()] = 0;

				output.push_back (Block (inner_block->source () + outer_block->source (), curr_source_idx[inner_block->source ()], curr_dest_idx, curr_size));

				commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION)
					<< "New block created: " << output.back () << std::endl
					<< "Inner source-block: " << *inner_block << std::endl
					<< "Outer source-block: " << *outer_block << std::endl;

				curr_source_idx[inner_block->source ()] += curr_size;
				curr_dest_idx += curr_size;
				rest_size = inner_block->size () - curr_size;

				if (rest_size == 0)
					++inner_block;
			}
		}
	}

	commentator.stop (MSG_DONE);

	return output;
}

Splicer &Splicer::compose (Splicer &output, const Splicer &inner, unsigned int horiz_inner_source, unsigned int vert_inner_source) const
{
	output._horiz_blocks.clear ();
	output._vert_blocks.clear ();

	map_blocks (output._horiz_blocks, _horiz_blocks, inner._horiz_blocks, horiz_inner_source);
	map_blocks (output._vert_blocks, _vert_blocks, inner._vert_blocks, vert_inner_source);

	return output;
}

bool Splicer::check_blocks (const std::vector<Block> blocks, const char *type) const
{
	std::ostringstream str;
	str << "Checking " << type << " blocks" << std::ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	std::vector<Block>::const_iterator i;
	size_t dest_idx = 0;
	std::map<unsigned int, size_t> src_idx;

	bool pass = true;

	for (i = blocks.begin (); i != blocks.end (); ++i) {
		if (dest_idx != i->destIndex ()) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Problem at block " << *i << ": Expected destination-index " << dest_idx << " but got " << i->destIndex () << std::endl;
			pass = false;
		}

		dest_idx = i->destIndex () + i->size ();

		if (src_idx.find (i->source ()) == src_idx.end () && i->sourceIndex () != 0) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Problem at block " << *i << ": Expected source-index 0 but got " << i->sourceIndex () << std::endl;
			pass = false;
		}
		else if (src_idx[i->source ()] != i->sourceIndex ()) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Problem at block " << *i << ": Expected source-index " << src_idx[i->source ()] << " but got " << i->sourceIndex () << std::endl;
			pass = false;
		}

		src_idx[i->source ()] = i->sourceIndex () + i->size ();
	}

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

bool Splicer::check () const
{
	bool pass = true;

	commentator.start ("Checking integrity of blocks in splicer", __FUNCTION__);

	pass = check_blocks (_horiz_blocks, "horizontal") && pass;
	pass = check_blocks (_vert_blocks, "vertical") && pass;

	commentator.stop (MSG_STATUS (pass));

	return pass;
}

std::ostream &operator << (std::ostream &os, const Splicer &splicer)
{
	std::vector<Block>::const_iterator i;

	os << "Horizontal blocks:" << std::endl;

	for (i = splicer._horiz_blocks.begin (); i != splicer._horiz_blocks.end (); ++i)
		os << "  " << *i << std::endl;

	os << "Vertical blocks:" << std::endl;

	for (i = splicer._vert_blocks.begin (); i != splicer._vert_blocks.end (); ++i)
		os << "  " << *i << std::endl;

	return os;
}

} // namespace LinBox

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
