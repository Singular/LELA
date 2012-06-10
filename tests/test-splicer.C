/* tests/test-splicer.C
 * Copyright 2011 Florian Fischer
 *
 * Written by Florian Fischer <florian-fischer@gmx.net>
 *
 * ---------------------------------------------------------
 *
 * See COPYING for license information.
 *
 * Generic test suite for splicer routines
 */

#include <iostream>

#include "lela/util/splicer.h"
#include "lela/util/commentator.h"
#include "lela/util/debug.h"
#include "lela/blas/context.h"
#include "lela/blas/level1.h"
#include "lela/vector/bit-vector.h"
#include "lela/vector/bit-subvector.h"
#include "lela/vector/sparse-subvector.h"
#include "lela/matrix/submatrix.h"

#include "test-common.h"

using namespace LELA;
using namespace std;

class grid1
{
public:
	typedef GridTypeNormal GridType;

	Block **horiz_array;
	Block **vert_array;
	size_t i, j, m, n;
	
	grid1(size_t _m, size_t _n)
	{
		m = _m;
		n = _n;
		horiz_array = new Block *[m];
		vert_array = new Block *[m];
		
		for(size_t k = 0; k < _m ; ++k)
		{
			horiz_array[k] = new Block [n];
			vert_array[k] = new Block [n];
		}
		i = j = 0;
	}
	
	void operator () (const Block &horiz_block, const Block &vert_block)
	{
		horiz_array[i][j] = horiz_block;
		vert_array[i][j] = vert_block;
		++j;
		
		if(j == n)
		{
			j = 0;
			++i;
		}
	}

};

class grid2
{
        bool &pass;
public:
	typedef GridTypeNormal GridType;

        grid2(bool &_pass):pass(_pass)
        {}

        void operator () (const Block &horiz_block, const Block &vert_block)
        {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: There are horizontal/vertical blocks, but they shoudn't." << endl;
		pass = false;
	}
};

bool testSplice ()
{
	bool pass = true;
        Block lasthoriBlock;
	Block lastvertiBlock;
   	grid1 g(4,4);

	std::ostringstream str;

	str << "Testing Splicer::splice " << std::ends;
	commentator.start (str.str ().c_str (), __FUNCTION__);

	Splicer S;
	S.addHorizontalBlock (Block (0, 0, 0, 0, 4));
	S.addHorizontalBlock (Block (0, 1, 4, 0, 4));
	S.addHorizontalBlock (Block (1, 0, 0, 4, 4));
	S.addHorizontalBlock (Block (1, 1, 4, 4, 4));
	S.addVerticalBlock (Block (0, 0, 0, 0, 4));
	S.addVerticalBlock (Block (0, 1, 4, 0, 4));
	S.addVerticalBlock (Block (1, 0, 0, 4, 4));
	S.addVerticalBlock (Block (1, 1, 4, 4, 4));

	S.splice(g);

	ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	for(size_t i = 0; i < 4; ++i)
	{
		for(size_t j = 0; j < 4; ++j)
                {
			report << "HorizontalBlock: " << g.horiz_array[i][j] << std::endl;
			report << "VerticalBlock: " << g.vert_array[i][j] << std::endl;
                
			if (g.horiz_array[i][j].dest () >= 2 )
			{
				report << "Destination id of the Horizontal Block: " << g.horiz_array[i][j].dest () << std::endl;
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR 1: Destiantion id out of bounds" << endl;
				pass = false;
			}

                        if (g.vert_array[i][j].dest () >= 2)
			{
				report << "Destination id of the Vertical Block: " << g.vert_array[i][j].dest () << std::endl;
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR 2: Destiantion id out of bounds" << endl;
				pass = false;
			}

			if (g.horiz_array[i][j].source () >= 2 )
			{
				report << "Source id of the Horizontal Block: " << g.horiz_array[i][j].source () << std::endl;
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR 3: Source id out of bounds" << endl;
				pass = false;
			}

                        if (g.vert_array[i][j].source () >= 2)
                        {
                                report << "Source id of the Vertical Block: " << g.vert_array[i][j].source () << std::endl;
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR 4: Source id out of bounds" << endl;
				pass = false;
			}

			if (g.horiz_array[i][j].sourceIndex () != 0 && g.horiz_array[i][j].sourceIndex () != 4)
			{
				report << "Source Index of the Horizontal Block: " << g.horiz_array[i][j].sourceIndex () << std::endl;
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR 5: Blocks do not allign" << endl;
				pass = false;
			}

			if ( g.vert_array[i][j].sourceIndex () != 0 && g.vert_array[i][j].sourceIndex () != 4)
			{
				report << "Source Index of the Vertical Block: " << g.vert_array[i][j].sourceIndex () << std::endl;
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR 6: Blocks do not allign" << endl;
				pass = false;
			}

			if (g.horiz_array[i][j].destIndex () != 0 && g.horiz_array[i][j].destIndex () != 4)
			{
				report << "Destination Index of the Horizontal Block: " << g.horiz_array[i][j].destIndex () << std::endl;
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR 7: Blocks do not allign" << endl;
				pass = false;
			}

			if ( g.vert_array[i][j].destIndex () != 0 && g.vert_array[i][j].destIndex () != 4 )
			{
				report << "Destination Index of the Vertical Block: " << g.vert_array[i][j].destIndex () << std::endl;
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR 8: Blocks do not allign" << endl;
				pass = false;
			}

			if (g.horiz_array[i][j].size () != 4 )
			{
				report << "Size of the Horizontal Block: " << g.vert_array[i][j].size () << std::endl;
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR 9: Incorrect size" << endl;
				pass = false;
			}

                        if (g.vert_array[i][j].size () != 4)
			{
                                report << "Size of the Vertical Block: " << g.vert_array[i][j].size () << std::endl;
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR 10: Incorrect size" << endl;
				pass = false;
			}

			if (g.horiz_array[i][j].dest () == 1 && g.horiz_array[i][j].sourceIndex () != 4)
			{
				report << "Destination id of the Horizontal Block" << g.horiz_array[i][j].dest () << std::endl;
				report << "Source Index of the Horizontal Block" << g.horiz_array[i][j].sourceIndex () << std::endl;
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR 11: Blocks do not allign" << endl;
				pass = false;
			}

                        if (g.vert_array[i][j].dest () == 1 && g.vert_array[i][j].sourceIndex () != 4)
                        {
                                report << "Destination id of the Vertical Block" << g.vert_array[i][j].dest () << std::endl;
                                report << "Source Index of the Vertical Block" << g.vert_array[i][j].sourceIndex () << std::endl;
                                commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR 12: Blocks do not allign" << endl;
				pass = false;
                        }

			if (g.horiz_array[i][j].source () == 0 && g.horiz_array[i][j].destIndex () != 0 )
			{
                                report << "Source id of the Horizontal Block" << g.horiz_array[i][j].source () << std::endl;
                                report << "Destination Index of the Horizontal Block" << g.horiz_array[i][j].destIndex () << std::endl;
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR 13: Blocks do not allign" << endl;
				pass = false;
			}

			if (g.vert_array[i][j].source () == 0 && g.vert_array[i][j].destIndex () != 0)
                        {
				report << "Source id of the Vertical Block" << g.vert_array[i][j].source () << std::endl;
                                report << "Destination Index of the Vertical Block" << g.vert_array[i][j].destIndex () << std::endl;
                                commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR 14: Blocks do not allign" << endl;
                                pass = false;
                        }
		 
			if(i == 0)
			{
				if(g.horiz_array[i][j].dest () != 0)
				{
					commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
						<< "ERROR 15: Blocks are in the wrong order" << endl;
					report << "Destination id of the Horizontal Block: " << g.horiz_array[i][j].dest () << std::endl;
					pass = false;
				}

				if(g.horiz_array[i][j].source () != 0)
				{
					commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
						<< "ERROR 16: Blocks are in the wrong order" << endl;
					report << "Source id of the Horizontal Block: " << g.horiz_array[i][j].source () << std::endl;
					pass = false;
				}

				if(g.horiz_array[i][j].destIndex () != 0)
				{
					commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
						<< "ERROR 17: Blocks are in the wrong order" << endl;
					report << "Destination Index of the Horizontal Block: " << g.horiz_array[i][j].destIndex () << std::endl;
					pass = false;
				}

				if(g.horiz_array[i][j].sourceIndex () != 0)
                                {
                                        commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
						<< "ERROR 18: Blocks are in the wrong order" << endl;
                                        report << "Source Index of the Horizontal Block: " << g.horiz_array[i][j].sourceIndex () << std::endl;
                                        pass = false;
                                }

				if(g.horiz_array[i][j].size () != 4)
                                {
                                        commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                                                << "ERROR 19: Blocks are in the wrong order" << endl;
                                        report << "Size of the Horizontal Block: " << g.horiz_array[i][j].size () << std::endl;
                                        pass = false;
                                }				

			}

			else
			{
				if ( (g.horiz_array[i-1][j].source () == 0 && g.horiz_array[i-1][j].dest () == 0) && 
				     (g.horiz_array[i][j].source () != 0 || g.horiz_array[i][j].dest () != 1) )
				{
					report << "Source id of the Previous Horizontal Block: " << g.horiz_array[i-1][j].source () << std::endl;
					report << "Destination id of the Previous Horizontal Block: " << g.horiz_array[i-1][j].dest () << std::endl;
                                        report << "Source id of the Horizontal Block: " << g.horiz_array[i][j].source () << std::endl;
                                        report << "Destination id of the Horizontal Block: " << g.horiz_array[i][j].dest () << std::endl;
					commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
						<< "ERROR 25: Blocks are in the wrong order" << endl;
					pass = false;
				}

				if ( (g.horiz_array[i-1][j].source () == 0 && g.horiz_array[i-1][j].dest () == 1) && 
				     (g.horiz_array[i][j].source () != 1 || g.horiz_array[i][j].dest () != 0) )
				{
					report << "Source id of the Previous Horizontal Block: " << g.horiz_array[i-1][j].source () << std::endl;
                                        report << "Destination id of the Previous Horizontal Block: " << g.horiz_array[i-1][j].dest () << std::endl;
                                        report << "Source id of the Horizontal Block: " << g.horiz_array[i][j].source () << std::endl;
                                        report << "Destination id of the Horizontal Block: " << g.horiz_array[i][j].dest () << std::endl;
					commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
						<< "ERROR 27: Blocks are in the wrong order" << endl;
					pass = false;
				}

				if ( (g.horiz_array[i-1][j].source () == 1 && g.horiz_array[i-1][j].dest () == 0) &&
				     (g.horiz_array[i][j].source () != 1 || g.horiz_array[i][j].dest () != 1) )
                                {
                                        report << "Source id of the Previous Horizontal Block: " << g.horiz_array[i-1][j].source () << std::endl;
                                        report << "Destination id of the Previous Horizontal Block: " << g.horiz_array[i-1][j].dest () << std::endl;
                                        report << "Source id of the Horizontal Block: " << g.horiz_array[i][j].source () << std::endl;
                                        report << "Destination id of the Horizontal Block: " << g.horiz_array[i][j].dest () << std::endl;
                                        commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                                                << "ERROR 29: Blocks are in the wrong order" << endl;
                                        pass = false;
                                }

			}

			if(j == 0)
			{
				if(g.vert_array[i][j].dest () != 0)
				{
					commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
						<< "ERROR 15: Blocks are in the wrong order" << endl;
					report << "Destination id of the Vertical Block: " << g.vert_array[i][j].dest () << std::endl;
					pass = false;
				}

				if(g.vert_array[i][j].source () != 0)
				{
					commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
						<< "ERROR 16: Blocks are in the wrong order" << endl;
					report << "Source id of the Vertical Block: " << g.vert_array[i][j].source () << std::endl;
					pass = false;
				}

				if(g.vert_array[i][j].destIndex () != 0)
				{
					commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
						<< "ERROR 17: Blocks are in the wrong order" << endl;
					report << "Destination Index of the Vertical Block: " << g.vert_array[i][j].destIndex () << std::endl;
					pass = false;
				}

				if(g.vert_array[i][j].sourceIndex () != 0)
                                {
                                        commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
						<< "ERROR 18: Blocks are in the wrong order" << endl;
                                        report << "Source Index of the Vertical Block: " << g.vert_array[i][j].sourceIndex () << std::endl;
                                        pass = false;
                                }

				if(g.vert_array[i][j].size () != 4)
                                {
                                        commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                                                << "ERROR 19: Blocks are in the wrong order" << endl;
                                        report << "Size of the Vertical Block: " << g.vert_array[i][j].size () << std::endl;
                                        pass = false;
                                }				

			}
			else
			{
				if ( (g.vert_array[i][j-1].source () == 0 && g.vert_array[i][j-1].dest () == 0) && 
				     (g.vert_array[i][j].source () != 0 || g.vert_array[i][j].dest () != 1) )
				{
					report << "Source id of the Previous Vertical Block: " << g.vert_array[i][j-1].source () << std::endl;
					report << "Destination id of the Previous Vertical Block: " << g.vert_array[i][j-1].dest () << std::endl;
                                        report << "Source id of the Vertical Block: " << g.vert_array[i][j].source () << std::endl;
                                        report << "Destination id of the Vertical Block: " << g.vert_array[i][j].dest () << std::endl;
					commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
						<< "ERROR 25: Blocks are in the wrong order" << endl;
					pass = false;
				}

				if ( (g.vert_array[i][j-1].source () == 0 && g.vert_array[i][j-1].dest () == 1) && 
				     (g.vert_array[i][j].source () != 1 || g.vert_array[i][j].dest () != 0) )
				{
					report << "Source id of the Previous Vertical Block: " << g.vert_array[i][j-1].source () << std::endl;
                                        report << "Destination id of the Previous Vertical Block: " << g.vert_array[i][j-1].dest () << std::endl;
                                        report << "Source id of the Vertical Block: " << g.vert_array[i][j].source () << std::endl;
                                        report << "Destination id of the Vertical Block: " << g.vert_array[i][j].dest () << std::endl;
					commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
						<< "ERROR 27: Blocks are in the wrong order" << endl;
					pass = false;
				}

				if ( (g.vert_array[1][j-1].source () == 1 && g.vert_array[i][j-1].dest () == 0) &&
				     (g.vert_array[i][j].source () != 1 || g.vert_array[i][j].dest () != 1) )
                                {
                                        report << "Source id of the Previous Vertical Block: " << g.vert_array[i][j-1].source () << std::endl;
                                        report << "Destination id of the Previous Vertical Block: " << g.vert_array[i][j-1].dest () << std::endl;
                                        report << "Source id of the Vertical Block: " << g.vert_array[i][j].source () << std::endl;
                                        report << "Destination id of the Vertical Block: " << g.vert_array[i][j].dest () << std::endl;
                                        commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                                                << "ERROR 29: Blocks are in the wrong order" << endl;
                                        pass = false;
                                }

			}


		}

	}

	commentator.stop (MSG_STATUS (pass));
	return pass;
}

bool testclearHorizontalBlocks ()
{
        bool pass = true;
	grid2 g (pass);

	std::ostringstream str;

        str << "Testing Splicer::clearHorizontalBlocks " << std::ends;
        commentator.start (str.str ().c_str (), __FUNCTION__);

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

	S.splice(g);

        commentator.stop (MSG_STATUS (pass));
        return pass;
}

bool testclearVerticalBlocks ()
{
        bool pass = true;
        grid2 g (pass);

	std::ostringstream str;

        str << "Testing Splicer::clearVerticalBlocks " << std::ends;
        commentator.start (str.str ().c_str (), __FUNCTION__);

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

        S.splice(g);

        commentator.stop (MSG_STATUS (pass));
        return pass;
}

bool testfillHorizontal ()
{
        bool pass = true;
        grid1 g (5,4);

	std::ostringstream str;

        str << "Testing Splicer::fillHorizontal " << std::ends;
        commentator.start (str.str ().c_str (), __FUNCTION__);

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

        S.splice(g);

        for(size_t i = 0; i < 4; ++i)
	{
		for(size_t j = 0; j < 4; ++j)
		{
				if (g.horiz_array[i][j].source () == 2)
				{				
					commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				       		<< "ERROR 2: Wrong block inserted." << endl;
			       		pass = false;
				}

				if (g.vert_array[i][j].source () == 2)
				{
					commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
						<< "ERROR 2.5: Wrong block inserted." << endl;
					pass = false;
				}

			}
			
		}

	for(size_t j = 0; j < 4; ++j)
	{
		if (g.horiz_array[4][j].source () != 2 || g.horiz_array[4][j].dest () != 0 || g.horiz_array[4][j].sourceIndex () !=  0 ||
		    g.horiz_array[4][j].destIndex () != 8 || g.horiz_array[4][j].size () != 12)
		{
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR 1: Wrong block inserted." << endl;
			pass = false;
		}

	}

	commentator.stop (MSG_STATUS (pass));
	return pass;
}

bool testfillVertical ()
{
        bool pass = true;
        grid1 g (4,5);

	std::ostringstream str;

        str << "Testing Splicer::fillVertical " << std::ends;
        commentator.start (str.str () .c_str (), __FUNCTION__);

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

        S.splice(g);

        for(size_t i = 0; i < 4; ++i)
	{
		for(size_t j = 0; j < 4; ++j)
		{
			if (g.horiz_array[i][j].source () == 2)
			{
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR 2: Wrong block inserted." << endl;
				pass = false;
			}

			if (g.vert_array[i][j].source () == 2)
			{
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR 2.5: Wrong block inserted." << endl;
				pass = false;
			}

		}

        }

        for(size_t i = 0; i < 4; ++i)
	{
		if (g.vert_array[i][4].source () != 2 || g.vert_array[i][4].dest () != 0 || g.vert_array[i][4].sourceIndex () !=  0 ||
		    g.vert_array[i][4].destIndex () != 8 || g.vert_array[i][4].size () != 12)
		{
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR 1: Wrong block inserted." << endl;
			pass = false;
		}

	}


        commentator.stop (MSG_STATUS (pass));
        return pass;
}

/*bool testConsolidate ()
{
       	bool pass = true;

	grid1 g (2,2);
	
	std::ostringstream str;

        str << "Testing Splicer::testconsolidate " << std::ends;
        commentator.start (str.str () .c_str (), __FUNCTION__);
	ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

        Splicer S;
        S.addHorizontalBlock (Block (0, 0, 0, 0, 2));
        S.addHorizontalBlock (Block (0, 0, 2, 0, 3));
        S.addHorizontalBlock (Block (0, 1, 0, 2, 2));
        S.addHorizontalBlock (Block (0, 1, 2, 2, 3));
        S.addVerticalBlock (Block (0, 0, 0, 0, 2));
        S.addVerticalBlock (Block (0, 0, 2, 0, 3));
        S.addVerticalBlock (Block (0, 1, 0, 2, 2));
        S.addVerticalBlock (Block (0, 1, 2, 2, 3));

        S.consolidate ();

        S.splice(g);

        for(size_t i = 0; i < 2; ++i)
	{
		for(size_t j = 0; j < 2; ++j)
		{
			if( i != 1)
			{
				if ( g.horiz_array[i][j].source () == g.horiz_array[i+1][j].source () )
				{
					if(g.horiz_array[i][j].dest () == g.horiz_array[i+1][j].dest () )
					{
						commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
							<< "ERROR 1: Consolidation of 2 horizontal blocks didn't happen." << endl;
						pass = false;
					}
				}
			}

			if( j != 1)
			{
				if (g.vert_array[i][j].source () == g.vert_array[i][j+1].source ())
				{
					if(g.vert_array[i][j].dest () == g.vert_array[i][j+1].dest () )
					{
						commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
							<< "ERROR 2: Consolidation of 2 vertical blocks didn't happen." << endl;
						pass = false;
					}
				}
			}
			
			if(g.horiz_array[i][j].source () != 0)
			{
				if(g.horiz_array[i][j].dest () != 0)
				{
					if(g.horiz_array[i][j].sourceIndex () != 0)
					{
						if(g.horiz_array[i][j].destIndex () != 0)
						{
							if(g.horiz_array[i][j].size () != 5)
							{
								commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
									<< "ERROR 3: Consolidation of 2 horizontal blocks went wrong." << endl;
								pass = false;
							}
						}
					}
				}
			}

			if(g.vert_array[i][j].source () != 0)
			{
				if(g.vert_array[i][j].dest () != 0)
				{
					if(g.vert_array[i][j].sourceIndex () != 0)
					{
						if(g.vert_array[i][j].destIndex () != 0)
						{
							if(g.vert_array[i][j].size () != 5)
							{
								commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
									<< "ERROR 4: Consolidation of 2 vertical blocks went wrong." << endl;
								pass = false;
							}
						}
					}
				}
			}
		}
        }

        commentator.stop (MSG_STATUS (pass));
        return pass;
	}*/

int main (int argc, char **argv)
{
	bool pass = true;

	static Argument args[] = {
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	commentator.start ("Testing Splicer", "Splicer");

	pass = testSplice () && pass;
	pass = testclearHorizontalBlocks () && pass;
	pass = testclearVerticalBlocks () && pass;
	pass = testfillHorizontal () && pass;
	pass = testfillVertical () && pass;
	// pass = testConsolidate () && pass;

	commentator.stop (MSG_STATUS (pass));
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: t
// c-basic-offset: 8
// End:

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
