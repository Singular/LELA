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
	linbox_check (!inner_blocks.empty ());

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
				output.push_back (Block (outer_block->source () + max_inner_block, outer_block->sourceIndex (), outer_block->destIndex (), outer_block->size ()));

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

class IsSourceP {
	unsigned int _source;

public:
	IsSourceP (int source) : _source (source) {}

	bool operator () (const Block &b)
		{ return b.source () == _source; }
};

void Splicer::fillHorizontal (unsigned int id, size_t rowdim)
{
	linbox_check (!_horiz_blocks.empty ());
	linbox_check (_horiz_blocks.back ().destIndex () + _horiz_blocks.back ().size () <= rowdim);

	if (_horiz_blocks.back ().destIndex () + _horiz_blocks.back ().size () == rowdim)
		return;

	std::vector<Block>::reverse_iterator i = std::find_if (_horiz_blocks.rbegin (), _horiz_blocks.rend (), IsSourceP (id));

	size_t source_idx;

	if (i == _horiz_blocks.rend ())
		source_idx = 0;
	else
		source_idx = i->sourceIndex () + i->size ();

	_horiz_blocks.push_back (Block (id, source_idx, _horiz_blocks.back ().destIndex () + _horiz_blocks.back ().size (),
					rowdim - (_horiz_blocks.back ().destIndex () + _horiz_blocks.back ().size ())));
}

void Splicer::fillVertical (unsigned int id, size_t coldim)
{
	linbox_check (!_vert_blocks.empty ());
	linbox_check (_vert_blocks.back ().destIndex () + _vert_blocks.back ().size () <= coldim);

	if (_vert_blocks.back ().destIndex () + _vert_blocks.back ().size () == coldim)
		return;

	std::vector<Block>::reverse_iterator i = std::find_if (_vert_blocks.rbegin (), _vert_blocks.rend (), IsSourceP (id));

	size_t source_idx;

	if (i == _vert_blocks.rend ())
		source_idx = 0;
	else
		source_idx = i->sourceIndex () + i->size ();

	_vert_blocks.push_back (Block (id, source_idx, _vert_blocks.back ().destIndex () + _vert_blocks.back ().size (),
				       coldim - (_vert_blocks.back ().destIndex () + _vert_blocks.back ().size ())));
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

void Splicer::substitute (std::vector<Block> &blocks, const std::vector<Block> &other_blocks, unsigned int source, unsigned int other_source)
{
	std::vector<Block>::iterator i;
	std::vector<Block>::const_iterator j = other_blocks.begin (), j_stop;
	size_t curr_source_idx = 0, curr_dest_idx = 0, size_from_source = 0;

	std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	for (i = blocks.begin (); i != blocks.end (); ++i) {
		if (i->source () == source) {
			for (j_stop = j; j_stop != other_blocks.end () && j_stop->destIndex () < i->sourceIndex () + i->size (); ++j_stop)
				if (j_stop->source () == other_source)
					size_from_source += j_stop->size ();

			if (size_from_source > 0) {
				report << "Replacing " << *i;
				size_t new_size = std::min (size_from_source, i->size ());
				*i = Block (source, curr_source_idx, curr_dest_idx, new_size);
				report << " with " << *i << " (substitution)" << std::endl;
				curr_dest_idx += new_size;
				curr_source_idx += new_size;
				size_from_source -= new_size;
			} else {
				report << "Erasing block " << *i << std::endl;
				i = blocks.erase (i);
				--i;
			}

			j = j_stop;
		} else {
			report << "Replacing " << *i;
			*i = Block (i->source (), i->sourceIndex (), curr_dest_idx, i->size ());
			report << " with " << *i << std::endl;
			curr_dest_idx += i->size ();
		}
	}
}

Splicer &Splicer::substituteHoriz (const Splicer &splicer, unsigned int source, unsigned int splicer_source)
{
	std::ostringstream str;
	str << "Substituting horizontal block " << splicer_source << " for block " << source << std::ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	report << "Splicer before substitution:" << std::endl << *this;
	report << "Splicer to be substituted:" << std::endl << splicer;

	substitute (_horiz_blocks, splicer._horiz_blocks, source, splicer_source);

	commentator.stop (MSG_DONE);

	return *this;
}

Splicer &Splicer::substituteVert (const Splicer &splicer, unsigned int source, unsigned int splicer_source)
{
	std::ostringstream str;
	str << "Substituting vertical block " << splicer_source << " for block " << source << std::ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	std::ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

	report << "Splicer before substitution:" << std::endl << *this;
	report << "Splicer to be substituted:" << std::endl << splicer;

	substitute (_vert_blocks, splicer._vert_blocks, source, splicer_source);

	commentator.stop (MSG_DONE);

	return *this;
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
