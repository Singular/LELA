/* lela/util/splicer.C
 * Copyright 2011 Bradford Hovinen <hovinen@gmail.com>
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * Class to chop matrices into pieces and splice them together
 *
 * ------------------------------------
 * 
 * This file is part of LELA, licensed under the GNU General Public
 * License version 3. See COPYING for more information.
 */

#include <iostream>
#include <sstream>
#include <algorithm>

#include "lela/util/splicer.h"
#include "lela/util/commentator.h"

namespace LELA
{

std::ostream &operator << (std::ostream &os, const Block &b)
{
	os << "Block (source=" << b.source () << ",dest=" << b.dest () << ",source-index=" << b.sourceIndex () << ",dest-index=" << b.destIndex () << ",size=" << b.size () << ")";

	return os;
}

struct CompareBlockSources { inline bool operator () (const Block &b1, const Block &b2) { return b1.source () < b2.source (); } };

std::vector<Block> &Splicer::map_blocks (std::vector<Block> &output,
					 const std::vector<Block> &outer_blocks,
					 const std::vector<Block> &inner_blocks,
					 unsigned int inner_source,
					 unsigned int only_source) const
{
	lela_check (!inner_blocks.empty ());

	std::vector<Block>::const_iterator outer_block;
	std::vector<Block>::const_iterator inner_block = inner_blocks.begin ();

	size_t rest_size = 0, curr_size = 0;
	std::map<size_t, size_t> curr_source_idx;
	std::map<size_t, size_t> curr_dest_idx;
	unsigned int use_source;

	commentator.start ("Mapping blocks", __FUNCTION__);

	size_t max_inner_block = (only_source == noSource) ? (std::max_element (inner_blocks.begin (), inner_blocks.end (), CompareBlockSources ()))->source () : 0;

	for (outer_block = outer_blocks.begin (); outer_block != outer_blocks.end (); ++outer_block) {
		if (curr_dest_idx.find (outer_block->dest ()) == curr_dest_idx.end ())
			curr_dest_idx[outer_block->dest ()] = 0;

		if (outer_block->source () != inner_source) {
			if (outer_block->source () < inner_source)
				output.push_back (*outer_block);
			else
				output.push_back (Block (outer_block->source () + max_inner_block, outer_block->dest (),
							 outer_block->sourceIndex (), outer_block->destIndex (), outer_block->size ()));

			curr_dest_idx[outer_block->dest ()] += outer_block->size ();

			commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION)
				<< "New block created (from outer block): " << output.back () << std::endl;
		} else {
			if (rest_size > 0) {
				curr_size = std::min (rest_size, outer_block->size ());

				if (only_source == noSource || inner_block->source () == only_source) {
					use_source = (only_source == noSource) ? (inner_block->source () + outer_block->source ()) : outer_block->source ();

					output.push_back (Block (use_source, outer_block->dest (), curr_source_idx[inner_block->source ()], curr_dest_idx[outer_block->dest ()], curr_size));

					commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION)
						<< "New block created (leftover): " << output.back () << std::endl
						<< "Inner source-block: " << *inner_block << std::endl
						<< "Outer source-block: " << *outer_block << std::endl;
				}

				curr_source_idx[inner_block->source ()] += curr_size;
				curr_dest_idx[outer_block->dest ()] += curr_size;
				rest_size -= curr_size;

				if (rest_size == 0)
					++inner_block;
			}

			while (inner_block != inner_blocks.end () && rest_size == 0 && outer_block->destIndex () + outer_block->size () > curr_dest_idx[outer_block->dest ()]) {
				curr_size = std::min (inner_block->size (), outer_block->destIndex () + outer_block->size () - curr_dest_idx[outer_block->dest ()]);

				if (curr_source_idx.find (inner_block->source ()) == curr_source_idx.end ())
					curr_source_idx[inner_block->source ()] = 0;

				if (only_source == noSource || inner_block->source () == only_source) {
					use_source = (only_source == noSource) ? (inner_block->source () + outer_block->source ()) : outer_block->source ();

					output.push_back (Block (use_source, outer_block->dest (), curr_source_idx[inner_block->source ()], curr_dest_idx[outer_block->dest ()], curr_size));

					commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION)
						<< "New block created: " << output.back () << std::endl
						<< "Inner source-block: " << *inner_block << std::endl
						<< "Outer source-block: " << *outer_block << std::endl;
				}

				curr_source_idx[inner_block->source ()] += curr_size;
				curr_dest_idx[outer_block->dest ()] += curr_size;
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

class IsDestP {
	unsigned int _dest;

public:
	IsDestP (int dest) : _dest (dest) {}

	bool operator () (const Block &b)
		{ return b.dest () == _dest; }
};

void Splicer::fill_blocks (std::vector<Block> &blocks, unsigned int sid, unsigned int did, size_t dim)
{
	lela_check (!blocks.empty ());

	std::vector<Block>::reverse_iterator i_s = std::find_if (blocks.rbegin (), blocks.rend (), IsSourceP (sid));
	std::vector<Block>::reverse_iterator i_d = std::find_if (blocks.rbegin (), blocks.rend (), IsDestP (did));

	size_t source_idx, dest_idx;

	if (i_s == blocks.rend ())
		source_idx = 0;
	else
		source_idx = i_s->sourceIndex () + i_s->size ();

	if (i_d == blocks.rend ())
		dest_idx = 0;
	else
		dest_idx = i_d->destIndex () + i_d->size ();

	lela_check (dest_idx <= dim);

	if (dest_idx < dim)
		blocks.push_back (Block (sid, did, source_idx, dest_idx, dim - dest_idx));
}

void Splicer::consolidate_blocks (std::vector<Block> &blocks)
{
	std::vector<Block>::iterator i = blocks.begin (), i_end;
	size_t size;

	while (i != blocks.end ()) {
		size = 0;

		for (i_end = i; i_end != blocks.end () && i_end->source () == i->source () && i_end->dest () == i->dest (); ++i_end)
			size += i_end->size ();

		*i = Block (i->source (), i->dest (), i->sourceIndex (), i->destIndex (), size);

		++i;
		i = blocks.erase (i, i_end);
	}
}

void Splicer::consolidate ()
{
	consolidate_blocks (_horiz_blocks);
	consolidate_blocks (_vert_blocks);
}

void Splicer::remove_gaps_from_blocks (std::vector<Block> &blocks)
{
	std::vector<Block>::iterator i;
	std::map<unsigned int, size_t> curr_source_idx, curr_dest_idx;

	for (i = blocks.begin (); i != blocks.end (); ++i) {
		if (curr_source_idx.find (i->source ()) == curr_source_idx.end ())
			curr_source_idx[i->source ()] = 0;

		if (curr_dest_idx.find (i->dest ()) == curr_dest_idx.end ())
			curr_dest_idx[i->dest ()] = 0;

		*i = Block (i->source (), i->dest (), curr_source_idx[i->source ()], curr_dest_idx[i->dest ()], i->size ());
		curr_source_idx[i->source ()] += i->size ();
		curr_dest_idx[i->dest ()] += i->size ();
	}
}

void Splicer::removeGaps ()
{
	remove_gaps_from_blocks (_horiz_blocks);
	remove_gaps_from_blocks (_vert_blocks);
}

Splicer &Splicer::compose (Splicer &output,
			   const Splicer &inner,
			   unsigned int horiz_inner_source,
			   unsigned int vert_inner_source,
			   unsigned int horiz_only_source,
			   unsigned int vert_only_source) const
{
	output._horiz_blocks.clear ();
	output._vert_blocks.clear ();

	map_blocks (output._horiz_blocks, _horiz_blocks, inner._horiz_blocks, horiz_inner_source, horiz_only_source);
	map_blocks (output._vert_blocks, _vert_blocks, inner._vert_blocks, vert_inner_source, vert_only_source);

	return output;
}

void Splicer::reverse_blocks (std::vector<Block> &out_blocks, const std::vector<Block> &in_blocks) const
{
	std::vector<Block>::const_iterator i;

	for (i = in_blocks.begin (); i != in_blocks.end (); ++i)
		out_blocks.push_back (Block (i->dest (), i->source (), i->destIndex (), i->sourceIndex (), i->size ()));
}

Splicer &Splicer::reverse (Splicer &output) const
{
	reverse_blocks (output._horiz_blocks, _horiz_blocks);
	reverse_blocks (output._vert_blocks, _vert_blocks);
	return output;
}

bool Splicer::check_blocks (const std::vector<Block> &blocks, const char *type) const
{
	std::ostringstream str;
	str << "Checking " << type << " blocks" << std::ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	std::vector<Block>::const_iterator i;
	std::map<unsigned int, size_t> src_idx;
	std::map<unsigned int, size_t> dest_idx;

	bool pass = true;

	for (i = blocks.begin (); i != blocks.end (); ++i) {
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

		if (dest_idx.find (i->dest ()) == dest_idx.end () && i->destIndex () != 0) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Problem at block " << *i << ": Expected dest-index 0 but got " << i->destIndex () << std::endl;
			pass = false;
		}
		else if (dest_idx[i->dest ()] != i->destIndex ()) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Problem at block " << *i << ": Expected dest-index " << dest_idx[i->dest ()] << " but got " << i->destIndex () << std::endl;
			pass = false;
		}

		src_idx[i->source ()] = i->sourceIndex () + i->size ();
		dest_idx[i->dest ()] = i->destIndex () + i->size ();
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

} // namespace LELA

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
