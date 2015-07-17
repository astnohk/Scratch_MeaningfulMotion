#include "Scratch_MeaningfulMotion.h"
#include "BRIEF.h"



bool
BRIEF(BRIEF *Features, const PNM &Img)
{
	ERROR Error("BRIEF");
	BRIEF_Patch Patch;

	if (Features == nullptr) {
		Error.Value("Features");
		Error.PointerNull();
		goto ExitError;
	}
	BRIEF_MakePatch(&Patch);
	return true;
// Error
ExitError:
	return false;
}

bool
BRIEF_MakePatch(BRIEF_Patch *Patch)
{
	ERROR Error("BRIEF_MakePatch");

	return true;
// Error
ExitError:
	return false;
}

