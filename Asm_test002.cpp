// Asm_test01.cpp : Defines the entry point for the console application.
//
#include "pch.h"
#include <immintrin.h> // MS version of immintrin.h covers AVX, AVX2 and FMA3

unsigned char RefPlane[100 * 100];
unsigned char CurrBlock[16 * 16];

unsigned int nSrcPitch[3] = { 32,16,16 }; // for 4 blocks 8x8
unsigned int nRefPitch[3] = { 100,100,100 };

int penaltyNew = 50;

#define SAD_4blocks \
/*load ref*/ \
ymm8_r1 = _mm256_loadu_si256((__m256i*)(pucRef + nRefPitch[0] * (y + 0) + x)); \
ymm9_r2 = _mm256_loadu_si256((__m256i*)(pucRef + nRefPitch[0] * (y + 1) + x)); \
ymm10_r3 = _mm256_loadu_si256((__m256i*)(pucRef + nRefPitch[0] * (y + 2) + x)); \
ymm11_r4 = _mm256_loadu_si256((__m256i*)(pucRef + nRefPitch[0] * (y + 3) + x)); \
ymm12_r5 = _mm256_loadu_si256((__m256i*)(pucRef + nRefPitch[0] * (y + 4) + x)); \
ymm13_r6 = _mm256_loadu_si256((__m256i*)(pucRef + nRefPitch[0] * (y + 5) + x)); \
ymm14_r7 = _mm256_loadu_si256((__m256i*)(pucRef + nRefPitch[0] * (y + 6) + x)); \
ymm15_r8 = _mm256_loadu_si256((__m256i*)(pucRef + nRefPitch[0] * (y + 7) + x)); \
/* calc sads */ \
ymm8_r1 = _mm256_sad_epu8(ymm8_r1, ymm0_src_r1); \
ymm9_r2 = _mm256_sad_epu8(ymm9_r2, ymm1_src_r2); \
ymm10_r3 = _mm256_sad_epu8(ymm10_r3, ymm2_src_r3); \
ymm11_r4 = _mm256_sad_epu8(ymm11_r4, ymm3_src_r4); \
ymm12_r5 = _mm256_sad_epu8(ymm12_r5, ymm4_src_r5); \
ymm13_r6 = _mm256_sad_epu8(ymm13_r6, ymm5_src_r6); \
ymm14_r7 = _mm256_sad_epu8(ymm14_r7, ymm6_src_r7); \
ymm15_r8 = _mm256_sad_epu8(ymm15_r8, ymm7_src_r8); \
\
ymm8_r1 = _mm256_adds_epu16(ymm8_r1, ymm9_r2); \
ymm10_r3 = _mm256_adds_epu16(ymm10_r3, ymm11_r4); \
ymm12_r5 = _mm256_adds_epu16(ymm12_r5, ymm13_r6); \
ymm14_r7 = _mm256_adds_epu16(ymm14_r7, ymm15_r8); \
\
ymm8_r1 = _mm256_adds_epu16(ymm8_r1, ymm10_r3); \
ymm12_r5 = _mm256_adds_epu16(ymm12_r5, ymm14_r7); \
\
ymm8_r1 = _mm256_adds_epu16(ymm8_r1, ymm12_r5);

#define _mm_cmpge_epu16(a, b) \
        _mm_cmpeq_epi16(_mm_max_epu16(a, b), a)


int main()
{
	unsigned char RefPlane[100 * 100];
	unsigned char CurrBlock[16 * 16];

	unsigned char CurrBlocks4[8 * 4 * 8];

	for (int i = 0; i < 16; i++)
	{
		for (int j = 0; j < 16; j++)
		{
			CurrBlock[j * 16 + i] = i + j * 10;
		}
	}

	for (int i = 0; i < 100; i++) // 
	{
		for (int j = 0; j < 100; j++) //
		{
			RefPlane[i * 100 + j] = 10 * (i - 2) + (j - 2);
		}
		int idbr02 = 0;
	}


	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 32; j++)
		{
			CurrBlocks4[i * 32 + j] = j + i * 10;
		}
	}


	unsigned char *pucRef = &RefPlane[0]; // upper left corner
	unsigned char *pucCurr = &CurrBlocks4[0];

	// 4 blocks at once proc, pel=1
	__m256i ymm0_src_r1 = _mm256_loadu_si256((__m256i*)(pucCurr + nSrcPitch[0] * 0));
	__m256i	ymm1_src_r2 = _mm256_loadu_si256((__m256i*)(pucCurr + nSrcPitch[0] * 1));
	__m256i	ymm2_src_r3 = _mm256_loadu_si256((__m256i*)(pucCurr + nSrcPitch[0] * 2));
	__m256i	ymm3_src_r4 = _mm256_loadu_si256((__m256i*)(pucCurr + nSrcPitch[0] * 3));
	__m256i	ymm4_src_r5 = _mm256_loadu_si256((__m256i*)(pucCurr + nSrcPitch[0] * 4));
	__m256i	ymm5_src_r6 = _mm256_loadu_si256((__m256i*)(pucCurr + nSrcPitch[0] * 5));
	__m256i	ymm6_src_r7 = _mm256_loadu_si256((__m256i*)(pucCurr + nSrcPitch[0] * 6));
	__m256i	ymm7_src_r8 = _mm256_loadu_si256((__m256i*)(pucCurr + nSrcPitch[0] * 7));

	__m256i ymm8_r1, ymm9_r2, ymm10_r3, ymm11_r4, ymm12_r5, ymm13_r6, ymm14_r7, ymm15_r8;

	int x, y;
	__m256i part_sads1, part_sads2;
#ifdef _DEBUG
	part_sads1 = _mm256_setzero_si256();
	part_sads2 = _mm256_setzero_si256();
#endif

	// 1st 4sads
	y = 0; x = 0;
	SAD_4blocks
/*			// load ref
	ymm8_r1 = _mm256_loadu_si256((__m256i*)(pucRef + nRefPitch[0] * (y + 0) + x));
	ymm9_r2 = _mm256_loadu_si256((__m256i*)(pucRef + nRefPitch[0] * (y + 1) + x));
	ymm10_r3 = _mm256_loadu_si256((__m256i*)(pucRef + nRefPitch[0] * (y + 2) + x));
	ymm11_r4 = _mm256_loadu_si256((__m256i*)(pucRef + nRefPitch[0] * (y + 3) + x));
	ymm12_r5 = _mm256_loadu_si256((__m256i*)(pucRef + nRefPitch[0] * (y + 4) + x));
	ymm13_r6 = _mm256_loadu_si256((__m256i*)(pucRef + nRefPitch[0] * (y + 5) + x));
	ymm14_r7 = _mm256_loadu_si256((__m256i*)(pucRef + nRefPitch[0] * (y + 6) + x));
	ymm15_r8 = _mm256_loadu_si256((__m256i*)(pucRef + nRefPitch[0] * (y + 7) + x));

	// calc sads
	ymm8_r1 = _mm256_sad_epu8(ymm8_r1, ymm0_src_r1);
	ymm9_r2 = _mm256_sad_epu8(ymm9_r2, ymm1_src_r2);
	ymm10_r3 = _mm256_sad_epu8(ymm10_r3, ymm2_src_r3);
	ymm11_r4 = _mm256_sad_epu8(ymm11_r4, ymm3_src_r4);
	ymm12_r5 = _mm256_sad_epu8(ymm12_r5, ymm4_src_r5);
	ymm13_r6 = _mm256_sad_epu8(ymm13_r6, ymm5_src_r6);
	ymm14_r7 = _mm256_sad_epu8(ymm14_r7, ymm6_src_r7);
	ymm15_r8 = _mm256_sad_epu8(ymm15_r8, ymm7_src_r8);

	ymm8_r1 = _mm256_adds_epu16(ymm8_r1, ymm9_r2);
	ymm10_r3 = _mm256_adds_epu16(ymm10_r3, ymm11_r4);
	ymm12_r5 = _mm256_adds_epu16(ymm12_r5, ymm13_r6);
	ymm14_r7 = _mm256_adds_epu16(ymm14_r7, ymm15_r8);

	ymm8_r1 = _mm256_adds_epu16(ymm8_r1, ymm10_r3);
	ymm12_r5 = _mm256_adds_epu16(ymm12_r5, ymm14_r7);

	ymm8_r1 = _mm256_adds_epu16(ymm8_r1, ymm12_r5);
	*/
	part_sads1 = _mm256_blend_epi16(part_sads1, ymm8_r1, 17);
	part_sads1 = _mm256_slli_si256(part_sads1, 2);

	// 2nd 4sads
	x = 1; y = 0;
	SAD_4blocks

	part_sads1 = _mm256_blend_epi16(part_sads1, ymm8_r1, 17);
	part_sads1 = _mm256_slli_si256(part_sads1, 2);

	// 3rd 4sads
	x = 2; y = 0;
	SAD_4blocks

	part_sads1 = _mm256_blend_epi16(part_sads1, ymm8_r1, 17);
	part_sads1 = _mm256_slli_si256(part_sads1, 2);

	// 4th 4sads
	x = 0; y = 1;
	SAD_4blocks

	part_sads1 = _mm256_blend_epi16(part_sads1, ymm8_r1, 17); // part_sads1 ready 4x4

	// 5th 4sads
	x = 1; y = 1;
	SAD_4blocks

	part_sads2 = _mm256_blend_epi16(part_sads2, ymm8_r1, 17);
	part_sads2 = _mm256_slli_si256(part_sads2, 2);

	// 6th 4sads
	x = 2; y = 1;
	SAD_4blocks

	part_sads2 = _mm256_blend_epi16(part_sads2, ymm8_r1, 17);
	part_sads2 = _mm256_slli_si256(part_sads2, 2);

	// 7th 4sads
	x = 0; y = 2;
	SAD_4blocks

	part_sads2 = _mm256_blend_epi16(part_sads2, ymm8_r1, 17);
	part_sads2 = _mm256_slli_si256(part_sads2, 2);

	// 8th 4sads
	x = 1; y = 2;
	SAD_4blocks

	part_sads2 = _mm256_blend_epi16(part_sads2, ymm8_r1, 17); // part_sads2 ready 4x4

	// 9th 4sads
	x = 2; y = 2;
	SAD_4blocks

	// 9th 4 sads in ymm8_r1
/*	unsigned int ui9SAD1 = _mm256_extract_epi16(ymm8_r1, 0);
	unsigned int ui9SAD2 = _mm256_extract_epi16(ymm8_r1, 3);
	unsigned int ui9SAD3 = _mm256_extract_epi16(ymm8_r1, 7);
	unsigned int ui9SAD4 = _mm256_extract_epi16(ymm8_r1, 15);
	*/
	__m256i ymm_tmp1, ymm_tmp2;
	// 8 SADs of 1 block
	ymm_tmp1 = _mm256_slli_si256(part_sads1, 8);
	ymm_tmp2 = _mm256_blend_epi16(part_sads2, ymm_tmp1, 240);

	__m128i xmm0_sad1 = _mm256_castsi256_si128(ymm_tmp2);

	// 8 SADS of 2 block
	ymm_tmp1 = _mm256_srli_si256(part_sads2, 8);
	ymm_tmp2 = _mm256_blend_epi16(part_sads1, ymm_tmp1, 15);

	__m128i xmm1_sad2 = _mm256_castsi256_si128(ymm_tmp2);

	part_sads1 = _mm256_permute4x64_epi64(part_sads1, 14); // move high 128bits to low 128 bits
	part_sads2 = _mm256_permute4x64_epi64(part_sads2, 14); // move high 128bits to low 128 bits

	// 8 SADs of 3 block
	ymm_tmp1 = _mm256_slli_si256(part_sads1, 8);
	ymm_tmp2 = _mm256_blend_epi16(part_sads2, ymm_tmp1, 240);

	__m128i xmm2_sad3 = _mm256_castsi256_si128(ymm_tmp2);
	
	// 8 SADs of 4 block
	ymm_tmp1 = _mm256_srli_si256(part_sads2, 8);
	ymm_tmp2 = _mm256_blend_epi16(part_sads1, ymm_tmp1, 15);

	__m128i xmm3_sad4 = _mm256_castsi256_si128(ymm_tmp2);

	unsigned int uiSADRes1 = _mm_cvtsi128_si32(_mm_minpos_epu16(xmm0_sad1));
	unsigned int uiSADRes2 = _mm_cvtsi128_si32(_mm_minpos_epu16(xmm1_sad2));
	unsigned int uiSADRes3 = _mm_cvtsi128_si32(_mm_minpos_epu16(xmm2_sad3));
	unsigned int uiSADRes4 = _mm_cvtsi128_si32(_mm_minpos_epu16(xmm3_sad4));

	
	unsigned short minsad = 65535;
	int idx_min_sad = 0;
/* SIMD it too
	if ((unsigned short)uiSADRes1 < minsad)
	{
		minsad = (unsigned short)uiSADRes1;
		idx_min_sad = 7 - (uiSADRes1 >> 16);
	}
	*/
	__m128i xmm_sad_res1234 = _mm_cvtsi32_si128(uiSADRes1); // replace with shift + blend
	xmm_sad_res1234 = _mm_insert_epi32(xmm_sad_res1234, uiSADRes2, 1);
	xmm_sad_res1234 = _mm_insert_epi32(xmm_sad_res1234, uiSADRes3, 2);
	xmm_sad_res1234 = _mm_insert_epi32(xmm_sad_res1234, uiSADRes4, 3);

	__m128i xmm_minsad1234 = xmm_sad_res1234;//_mm_set1_epi16(-1);
	__m128i idx_minsad1234 = _mm_set1_epi16(7);

//	xmm_minsad = _mm_min_epu16(xmm_minsad, xmm_sad_res1234);
	idx_minsad1234 = _mm_subs_epu16(idx_minsad1234, xmm_sad_res1234);
/*
	if ((unsigned short)ui9SAD1 < minsad)
	{
		minsad = (unsigned short)ui9SAD1 < minsad;
		idx_min_sad = 8;
	}
	*/
/*	__m128i xmm_sad_res1234_9 = _mm_cvtsi32_si128(ui9SAD1); // repack ymm8_r1 !!!
	xmm_sad_res1234_9 = _mm_insert_epi32(xmm_sad_res1234, ui9SAD2, 1);
	xmm_sad_res1234_9 = _mm_insert_epi32(xmm_sad_res1234, ui9SAD3, 2);
	xmm_sad_res1234_9 = _mm_insert_epi32(xmm_sad_res1234, ui9SAD4, 3);
	*/
	__m128i xmm_sad_res1234_9 = _mm256_castsi256_si128(_mm256_permutevar8x32_epi32(ymm8_r1, _mm256_setr_epi32(7, 7, 7, 7, 6, 4, 2, 0))); // check

	__m128i idx_minsad8 = _mm_set1_epi16(8);
	__m128i mask_idx8 = _mm_cmpge_epu16(xmm_sad_res1234, xmm_sad_res1234_9);

	idx_minsad1234 = _mm_blendv_epi8(idx_minsad1234, idx_minsad8, mask_idx8);
	xmm_minsad1234 = _mm_min_epu16(xmm_minsad1234, xmm_sad_res1234_9);

	int idx_min_sad1 = _mm_extract_epi16(idx_minsad1234, 0);
	int idx_min_sad2 = _mm_extract_epi16(idx_minsad1234, 2);
	int idx_min_sad3 = _mm_extract_epi16(idx_minsad1234, 4);
	int idx_min_sad4 = _mm_extract_epi16(idx_minsad1234, 6);

	int x_minsad1 = (idx_min_sad1 % 3) - 1; //- just comment where from x,y minsad come from
	int y_minsad1 = (idx_min_sad1 / 3) - 1;

	int x_minsad2 = (idx_min_sad2 % 3) - 1; //- just comment where from x,y minsad come from
	int y_minsad2 = (idx_min_sad2 / 3) - 1;

	int x_minsad3 = (idx_min_sad3 % 3) - 1; //- just comment where from x,y minsad come from
	int y_minsad3 = (idx_min_sad3 / 3) - 1;

	int x_minsad4 = (idx_min_sad4 % 3) - 1; //- just comment where from x,y minsad come from
	int y_minsad4 = (idx_min_sad4 / 3) - 1;

	int minsad1 = _mm_extract_epi16(idx_minsad1234, 1);
	int minsad2 = _mm_extract_epi16(idx_minsad1234, 3);
	int minsad3 = _mm_extract_epi16(idx_minsad1234, 5);
	int minsad4 = _mm_extract_epi16(idx_minsad1234, 7);

	/*	
	sad_t cost = minsad + ((penaltyNew * minsad) >> 8);
	if (cost >= workarea.nMinCost)
	{
		_mm256_zeroupper();
		return;
	}*/

	__m128i xmm_penaltyNew = _mm_set1_epi32(penaltyNew);
	__m128i xmm_minsad32 = _mm_srli_epi32(xmm_minsad1234, 16);
	xmm_minsad32 = _mm_mul_epi32(xmm_minsad32, xmm_penaltyNew);
	xmm_minsad32 = _mm_srli_epi32(xmm_minsad32, 8);
	int cost1 = _mm_extract_epi16(idx_minsad1234, 1);
	int cost2 = _mm_extract_epi16(idx_minsad1234, 3);
	int cost3 = _mm_extract_epi16(idx_minsad1234, 5);
	int cost4 = _mm_extract_epi16(idx_minsad1234, 7);


/*
	workarea.bestMV.x = mvx + (idx_min_sad % 3) - 1;
	workarea.bestMV.y = mvy + (idx_min_sad / 3) - 1;
	workarea.nMinCost = cost;
	workarea.bestMV.sad = minsad;
	*/




	return 0;
}

