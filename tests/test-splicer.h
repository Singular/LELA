/* tests/test-splicer.h
 * Copyright 2011 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@gmail.com>
 *
 * ---------------------------------------------------------
 *
 * See COPYING for license information.
 *
 * Generic test suite for splicer routines
 */

#ifndef __LELA_UTIL_SPLICER_TCC
#define __LELA_UTIL_SPLICER_TCC

#include "lela/util/splicer.h"
#include "lela/util/commentator.h"
#include "lela/util/debug.h"
#include "lela/blas/context.h"
#include "lela/blas/level1.h"
#include "lela/vector/bit-vector.h"
#include "lela/vector/bit-subvector.h"
#include "lela/vector/sparse-subvector.h"
#include "lela/matrix/submatrix.h"

class grid1
{
	bool &pass;
	Block lasthoriBlock;
	Block lastvertiBlock;
	bool firsttime;

public:
	int counter;
	grid1(bool &_pass):pass(_pass),firsttime(true), counter(0);
	{}

	void operator () (const Block &horiz_block, const Block &vert_block)
	{
		ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

		report << "HorizontalBlock: " << block << std::endl;

		report << "VerticalBlock: " << block << std::endl;
                
		if (horiz_block.dest () >= 2 || vert_block.dest () >= 2)
		{
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
       				<< "ERROR: Destiantion adress out of bounds" << endl;
	       		pass = false;
		}

		if (horiz_block.source () >= 2 || vert_block.source () >= 2)
		{
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Source adress out of bounds" << endl;
			pass = false;
		}

		if (horiz_block.src_idx () != 0 || horiz_block.src_idx () != 4 || vert_block.src_idx () != 0 ||  vert_block.src_idx () != 4)
		{
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Overlapping (gapping) blocks" << endl;
			pass = false;
		}

                if (horiz_block.dest_idx () != 0 || horiz_block.dest_idx () != 4 || vert_block.dest_idx () != 0 ||  vert_block.dest_idx () != 4)
		{	
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Overlapping (gapping) blocks" << endl;
			pass = false;
		}

                if (horiz_block.size () != 4 || vert_block.size () != 4)
		{
                        commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                                << "ERROR: Incorrect size" << endl;
                        pass = false;
		}

		if ((horiz_block.dest () == 1 && horiz_block.src_idx () != 4) ||(vert_block.dest () == 1 && vert_block.src_idx () != 4))
		{
                        commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                                << "ERROR: Overlapping (gapping) blocks" << endl;
                        pass = false;
		}

		if ((horiz_block.source () == 0 && horiz_block.dest_idx () != 0) ||(vert_block.source () == 0 && vert_block.dest_idx () != 0))
		{
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                                << "ERROR: Overlapping (gapping) blocks" << endl;
			pass = false;
		}
		if (firsttime)
		{
			firsttime = false;
		}
		else
		{
			if ( (lasthoriBlock.source () == 0 && lasthoriBlock.dest () == 0) && (horiz_block.source () != 0 || horiz_block.dest () != 1) )
			{
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR: Blocks are in the wrong order" << endl;
				pass = false;
			}

			if ( (lastvertiBlock.source () == 0 && lastvertiBlock.dest () == 0) && (vert_block.source () != 0 || vert_block.dest () != 1) )
			{
                                commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                                        << "ERROR: Blocks are in the wrong order" << endl;
				pass = false;
			}

			if ( (lasthoriBlock.source () == 0 && lasthoriBlock.dest () == 1) && (horiz_block.source () != 1 || horiz_block.dest () != 0) )
			{
                                commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                                        << "ERROR: Blocks are in the wrong order" << endl;
				pass = false;
			}

			if ( (lastvertiBlock.source () == 0 && lastvertiBlock.dest () == 1) && (vert_block.source () != 1 || vert_block.dest () != 0) )
			{
                                commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                                        << "ERROR: Blocks are in the wrong order" << endl;
				pass = false;
			}

			if ( (lasthoriBlock.source () == 1 && lasthoriBlock.dest () == 0) && (horiz_block.source () != 1 || horiz_block.dest () != 1) )
			{
                                commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                                        << "ERROR: Blocks are in the wrong order" << endl;
				pass = false;
			}

			if ( (lastvertiBlock.source () == 1 && lastvertiBlock.dest () == 0) && (vert_block.source () != 1 || vert_block.dest () != 1) )
			{
                                commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                                        << "ERROR: Blocks are in the wrong order" << endl;
				pass = false;
			}
		}

		lasthoriBlock = horiz_block;
		lastvertiBlock = vert_block;
		counter++;
	}
};

class grid2
{
        bool &pass;
public:
        grid2(bool &_pass):pass(_pass);
        {}

        void operator () (const Block &horiz_block, const Block &vert_block)
        {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: There are horizontal/vertical blocks, but they shoudn't." << endl;
		pass = false;
	}
};

class grid3
{
        bool &pass;
public:
	int counter;
        grid3(bool &_pass):pass(_pass),counter(0);
        {}

        void operator () (const Block &horiz_block, const Block &vert_block)
        {
                if (counter >= 16)
		{
			if (horiz_block.source () != 2 || horiz_block.dest () != 0 || horiz_block.src_idx () !=  0 || horiz_block.dest_idx () != 8 || horiz_block.size () != 4)
			{
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR: Wrong block inserted." << endl;
				pass = false;
			}
	        }
		else
		{
			if (horiz_block.source () == 2)
				{
					commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
						<< "ERROR: Wrong block inserted." << endl;
					pass = false;
				}
		}
        }
};


bool testSplice ()
{
	bool pass = true;
	grid1 g;

	std::ostringstream str;

	str << "Testing " << text << " Splicer::splice " << std::ends;
	commentator.start (str.str ().c_str (), __FUNCTION__, stream1.size ());

	Splicer S;
	S.addHorizontalBlock (Block (0, 0, 0, 0, 4));
	S.addHorizontalBlock (Block (0, 1, 4, 0, 4));
	S.addHorizontalBlock (Block (1, 0, 0, 4, 4));
	S.addHorizontalBlock (Block (1, 1, 4, 4, 4));
	S.addVerticalBlock (Block (0, 0, 0, 0, 4));
	S.addVerticalBlock (Block (0, 1, 4, 0, 4));
	S.addVerticalBlock (Block (1, 0, 0, 4, 4));
	S.addVerticalBlock (Block (1, 1, 4, 4, 4));

	Splicer::splice(g);

	if (g.counter != 16) pass = false;

	commentator.stop (MSG_STATUS (pass));
	return pass;
}

bool testclearHorizontalBlocks ()
{
        bool pass = true;
	grid2 g;

	std::ostringstream str;

        str << "Testing " << text << " Splicer::clearHorizontalBlocks " << std::ends;
        commentator.start (str.str ().c_str (), __FUNCTION__, stream1.size ());

        Splicer S;
        S.addHorizontalBlock (Block (0, 0, 0, 0, 4));
        S.addHorizontalBlock (Block (0, 1, 4, 0, 4));
        S.addHorizontalBlock (Block (1, 0, 0, 4, 4));
        S.addHorizontalBlock (Block (1, 1, 4, 4, 4));
        S.addVerticalBlock (Block (0, 0, 0, 0, 4));
        S.addVerticalBlock (Block (0, 1, 4, 0, 4));
        S.addVerticalBlock (Block (1, 0, 0, 4, 4));
        S.addVerticalBlock (Block (1, 1, 4, 4, 4));

	S.clearHorizontalBlocks ();

	Splicer.splice(g);

        commentator.stop (MSG_STATUS (pass));
        return pass;
}

bool testclearVerticalBlocks ()
{
        bool pass = true;
        grid2 g;

	std::ostringstream str;

        str << "Testing " << text << " Splicer::clearVerticalBlocks " << std::ends;
        commentator.start (str.str ().c_str (), __FUNCTION__, stream1.size ());

        Splicer S;
        S.addHorizontalBlock (Block (0, 0, 0, 0, 4));
        S.addHorizontalBlock (Block (0, 1, 4, 0, 4));
        S.addHorizontalBlock (Block (1, 0, 0, 4, 4));
        S.addHorizontalBlock (Block (1, 1, 4, 4, 4));
        S.addVerticalBlock (Block (0, 0, 0, 0, 4));
        S.addVerticalBlock (Block (0, 1, 4, 0, 4));
        S.addVerticalBlock (Block (1, 0, 0, 4, 4));
        S.addVerticalBlock (Block (1, 1, 4, 4, 4));

        S.clearVerticalBlocks ();

        Splicer.splice(g);

        commentator.stop (MSG_STATUS (pass));
        return pass;
}

bool testfillHorizontal ()
{
        bool pass = true;
        grid3 g;

	std::ostringstream str;

        str << "Testing " << text << " Splicer::fillHorizontal " << std::ends;
        commentator.start (str.str ().c_str (), __FUNCTION__, stream1.size ());

        Splicer S;
        S.addHorizontalBlock (Block (0, 0, 0, 0, 4));
        S.addHorizontalBlock (Block (0, 1, 4, 0, 4));
        S.addHorizontalBlock (Block (1, 0, 0, 4, 4));
        S.addHorizontalBlock (Block (1, 1, 4, 4, 4));
        S.addVerticalBlock (Block (0, 0, 0, 0, 4));
        S.addVerticalBlock (Block (0, 1, 4, 0, 4));
        S.addVerticalBlock (Block (1, 0, 0, 4, 4));
        S.addVerticalBlock (Block (1, 1, 4, 4, 4));

        S.fillHorizontal (2, 0, 20);

        Splicer.splice(g);

        commentator.stop (MSG_STATUS (pass));
        return pass;
}

bool testfillVertical ()
{
        bool pass = true;
        grid3 g;

	std::ostringstream str;

        str << "Testing " << text << " Splicer::fillVertical " << std::ends;
        commentator.start (str.str ().c_str (), __FUNCTION__, stream1.size ());

        Splicer S;
        S.addHorizontalBlock (Block (0, 0, 0, 0, 4));
        S.addHorizontalBlock (Block (0, 1, 4, 0, 4));
        S.addHorizontalBlock (Block (1, 0, 0, 4, 4));
        S.addHorizontalBlock (Block (1, 1, 4, 4, 4));
        S.addVerticalBlock (Block (0, 0, 0, 0, 4));
        S.addVerticalBlock (Block (0, 1, 4, 0, 4));
        S.addVerticalBlock (Block (1, 0, 0, 4, 4));
        S.addVerticalBlock (Block (1, 1, 4, 4, 4));

        S.fillVertical (2, 0, 20);

        Splicer.splice(g);

        commentator.stop (MSG_STATUS (pass));
        return pass;
}


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
