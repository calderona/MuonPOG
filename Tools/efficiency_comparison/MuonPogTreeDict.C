// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME MuonPogTreeDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "../src/MuonPogTree.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_muon_pogcLcLGenInfo(void *p = 0);
   static void *newArray_muon_pogcLcLGenInfo(Long_t size, void *p);
   static void delete_muon_pogcLcLGenInfo(void *p);
   static void deleteArray_muon_pogcLcLGenInfo(void *p);
   static void destruct_muon_pogcLcLGenInfo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::muon_pog::GenInfo*)
   {
      ::muon_pog::GenInfo *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::muon_pog::GenInfo >(0);
      static ::ROOT::TGenericClassInfo 
         instance("muon_pog::GenInfo", ::muon_pog::GenInfo::Class_Version(), "../src/MuonPogTree.h", 12,
                  typeid(::muon_pog::GenInfo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::muon_pog::GenInfo::Dictionary, isa_proxy, 4,
                  sizeof(::muon_pog::GenInfo) );
      instance.SetNew(&new_muon_pogcLcLGenInfo);
      instance.SetNewArray(&newArray_muon_pogcLcLGenInfo);
      instance.SetDelete(&delete_muon_pogcLcLGenInfo);
      instance.SetDeleteArray(&deleteArray_muon_pogcLcLGenInfo);
      instance.SetDestructor(&destruct_muon_pogcLcLGenInfo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::muon_pog::GenInfo*)
   {
      return GenerateInitInstanceLocal((::muon_pog::GenInfo*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::muon_pog::GenInfo*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_muon_pogcLcLGenParticle(void *p = 0);
   static void *newArray_muon_pogcLcLGenParticle(Long_t size, void *p);
   static void delete_muon_pogcLcLGenParticle(void *p);
   static void deleteArray_muon_pogcLcLGenParticle(void *p);
   static void destruct_muon_pogcLcLGenParticle(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::muon_pog::GenParticle*)
   {
      ::muon_pog::GenParticle *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::muon_pog::GenParticle >(0);
      static ::ROOT::TGenericClassInfo 
         instance("muon_pog::GenParticle", ::muon_pog::GenParticle::Class_Version(), "../src/MuonPogTree.h", 23,
                  typeid(::muon_pog::GenParticle), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::muon_pog::GenParticle::Dictionary, isa_proxy, 4,
                  sizeof(::muon_pog::GenParticle) );
      instance.SetNew(&new_muon_pogcLcLGenParticle);
      instance.SetNewArray(&newArray_muon_pogcLcLGenParticle);
      instance.SetDelete(&delete_muon_pogcLcLGenParticle);
      instance.SetDeleteArray(&deleteArray_muon_pogcLcLGenParticle);
      instance.SetDestructor(&destruct_muon_pogcLcLGenParticle);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::muon_pog::GenParticle*)
   {
      return GenerateInitInstanceLocal((::muon_pog::GenParticle*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::muon_pog::GenParticle*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_muon_pogcLcLMETs(void *p = 0);
   static void *newArray_muon_pogcLcLMETs(Long_t size, void *p);
   static void delete_muon_pogcLcLMETs(void *p);
   static void deleteArray_muon_pogcLcLMETs(void *p);
   static void destruct_muon_pogcLcLMETs(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::muon_pog::METs*)
   {
      ::muon_pog::METs *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::muon_pog::METs >(0);
      static ::ROOT::TGenericClassInfo 
         instance("muon_pog::METs", ::muon_pog::METs::Class_Version(), "../src/MuonPogTree.h", 47,
                  typeid(::muon_pog::METs), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::muon_pog::METs::Dictionary, isa_proxy, 4,
                  sizeof(::muon_pog::METs) );
      instance.SetNew(&new_muon_pogcLcLMETs);
      instance.SetNewArray(&newArray_muon_pogcLcLMETs);
      instance.SetDelete(&delete_muon_pogcLcLMETs);
      instance.SetDeleteArray(&deleteArray_muon_pogcLcLMETs);
      instance.SetDestructor(&destruct_muon_pogcLcLMETs);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::muon_pog::METs*)
   {
      return GenerateInitInstanceLocal((::muon_pog::METs*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::muon_pog::METs*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_muon_pogcLcLChambMatch(void *p = 0);
   static void *newArray_muon_pogcLcLChambMatch(Long_t size, void *p);
   static void delete_muon_pogcLcLChambMatch(void *p);
   static void deleteArray_muon_pogcLcLChambMatch(void *p);
   static void destruct_muon_pogcLcLChambMatch(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::muon_pog::ChambMatch*)
   {
      ::muon_pog::ChambMatch *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::muon_pog::ChambMatch >(0);
      static ::ROOT::TGenericClassInfo 
         instance("muon_pog::ChambMatch", ::muon_pog::ChambMatch::Class_Version(), "../src/MuonPogTree.h", 61,
                  typeid(::muon_pog::ChambMatch), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::muon_pog::ChambMatch::Dictionary, isa_proxy, 4,
                  sizeof(::muon_pog::ChambMatch) );
      instance.SetNew(&new_muon_pogcLcLChambMatch);
      instance.SetNewArray(&newArray_muon_pogcLcLChambMatch);
      instance.SetDelete(&delete_muon_pogcLcLChambMatch);
      instance.SetDeleteArray(&deleteArray_muon_pogcLcLChambMatch);
      instance.SetDestructor(&destruct_muon_pogcLcLChambMatch);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::muon_pog::ChambMatch*)
   {
      return GenerateInitInstanceLocal((::muon_pog::ChambMatch*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::muon_pog::ChambMatch*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_muon_pogcLcLHitInfo(void *p = 0);
   static void *newArray_muon_pogcLcLHitInfo(Long_t size, void *p);
   static void delete_muon_pogcLcLHitInfo(void *p);
   static void deleteArray_muon_pogcLcLHitInfo(void *p);
   static void destruct_muon_pogcLcLHitInfo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::muon_pog::HitInfo*)
   {
      ::muon_pog::HitInfo *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::muon_pog::HitInfo >(0);
      static ::ROOT::TGenericClassInfo 
         instance("muon_pog::HitInfo", ::muon_pog::HitInfo::Class_Version(), "../src/MuonPogTree.h", 84,
                  typeid(::muon_pog::HitInfo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::muon_pog::HitInfo::Dictionary, isa_proxy, 4,
                  sizeof(::muon_pog::HitInfo) );
      instance.SetNew(&new_muon_pogcLcLHitInfo);
      instance.SetNewArray(&newArray_muon_pogcLcLHitInfo);
      instance.SetDelete(&delete_muon_pogcLcLHitInfo);
      instance.SetDeleteArray(&deleteArray_muon_pogcLcLHitInfo);
      instance.SetDestructor(&destruct_muon_pogcLcLHitInfo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::muon_pog::HitInfo*)
   {
      return GenerateInitInstanceLocal((::muon_pog::HitInfo*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::muon_pog::HitInfo*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_muon_pogcLcLMuonFit(void *p = 0);
   static void *newArray_muon_pogcLcLMuonFit(Long_t size, void *p);
   static void delete_muon_pogcLcLMuonFit(void *p);
   static void deleteArray_muon_pogcLcLMuonFit(void *p);
   static void destruct_muon_pogcLcLMuonFit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::muon_pog::MuonFit*)
   {
      ::muon_pog::MuonFit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::muon_pog::MuonFit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("muon_pog::MuonFit", ::muon_pog::MuonFit::Class_Version(), "../src/MuonPogTree.h", 104,
                  typeid(::muon_pog::MuonFit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::muon_pog::MuonFit::Dictionary, isa_proxy, 4,
                  sizeof(::muon_pog::MuonFit) );
      instance.SetNew(&new_muon_pogcLcLMuonFit);
      instance.SetNewArray(&newArray_muon_pogcLcLMuonFit);
      instance.SetDelete(&delete_muon_pogcLcLMuonFit);
      instance.SetDeleteArray(&deleteArray_muon_pogcLcLMuonFit);
      instance.SetDestructor(&destruct_muon_pogcLcLMuonFit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::muon_pog::MuonFit*)
   {
      return GenerateInitInstanceLocal((::muon_pog::MuonFit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::muon_pog::MuonFit*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_muon_pogcLcLMuon(void *p = 0);
   static void *newArray_muon_pogcLcLMuon(Long_t size, void *p);
   static void delete_muon_pogcLcLMuon(void *p);
   static void deleteArray_muon_pogcLcLMuon(void *p);
   static void destruct_muon_pogcLcLMuon(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::muon_pog::Muon*)
   {
      ::muon_pog::Muon *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::muon_pog::Muon >(0);
      static ::ROOT::TGenericClassInfo 
         instance("muon_pog::Muon", ::muon_pog::Muon::Class_Version(), "../src/MuonPogTree.h", 133,
                  typeid(::muon_pog::Muon), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::muon_pog::Muon::Dictionary, isa_proxy, 4,
                  sizeof(::muon_pog::Muon) );
      instance.SetNew(&new_muon_pogcLcLMuon);
      instance.SetNewArray(&newArray_muon_pogcLcLMuon);
      instance.SetDelete(&delete_muon_pogcLcLMuon);
      instance.SetDeleteArray(&deleteArray_muon_pogcLcLMuon);
      instance.SetDestructor(&destruct_muon_pogcLcLMuon);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::muon_pog::Muon*)
   {
      return GenerateInitInstanceLocal((::muon_pog::Muon*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::muon_pog::Muon*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_muon_pogcLcLHLTObject(void *p = 0);
   static void *newArray_muon_pogcLcLHLTObject(Long_t size, void *p);
   static void delete_muon_pogcLcLHLTObject(void *p);
   static void deleteArray_muon_pogcLcLHLTObject(void *p);
   static void destruct_muon_pogcLcLHLTObject(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::muon_pog::HLTObject*)
   {
      ::muon_pog::HLTObject *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::muon_pog::HLTObject >(0);
      static ::ROOT::TGenericClassInfo 
         instance("muon_pog::HLTObject", ::muon_pog::HLTObject::Class_Version(), "../src/MuonPogTree.h", 278,
                  typeid(::muon_pog::HLTObject), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::muon_pog::HLTObject::Dictionary, isa_proxy, 4,
                  sizeof(::muon_pog::HLTObject) );
      instance.SetNew(&new_muon_pogcLcLHLTObject);
      instance.SetNewArray(&newArray_muon_pogcLcLHLTObject);
      instance.SetDelete(&delete_muon_pogcLcLHLTObject);
      instance.SetDeleteArray(&deleteArray_muon_pogcLcLHLTObject);
      instance.SetDestructor(&destruct_muon_pogcLcLHLTObject);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::muon_pog::HLTObject*)
   {
      return GenerateInitInstanceLocal((::muon_pog::HLTObject*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::muon_pog::HLTObject*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_muon_pogcLcLL1Muon(void *p = 0);
   static void *newArray_muon_pogcLcLL1Muon(Long_t size, void *p);
   static void delete_muon_pogcLcLL1Muon(void *p);
   static void deleteArray_muon_pogcLcLL1Muon(void *p);
   static void destruct_muon_pogcLcLL1Muon(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::muon_pog::L1Muon*)
   {
      ::muon_pog::L1Muon *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::muon_pog::L1Muon >(0);
      static ::ROOT::TGenericClassInfo 
         instance("muon_pog::L1Muon", ::muon_pog::L1Muon::Class_Version(), "../src/MuonPogTree.h", 293,
                  typeid(::muon_pog::L1Muon), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::muon_pog::L1Muon::Dictionary, isa_proxy, 4,
                  sizeof(::muon_pog::L1Muon) );
      instance.SetNew(&new_muon_pogcLcLL1Muon);
      instance.SetNewArray(&newArray_muon_pogcLcLL1Muon);
      instance.SetDelete(&delete_muon_pogcLcLL1Muon);
      instance.SetDeleteArray(&deleteArray_muon_pogcLcLL1Muon);
      instance.SetDestructor(&destruct_muon_pogcLcLL1Muon);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::muon_pog::L1Muon*)
   {
      return GenerateInitInstanceLocal((::muon_pog::L1Muon*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::muon_pog::L1Muon*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_muon_pogcLcLHLT(void *p = 0);
   static void *newArray_muon_pogcLcLHLT(Long_t size, void *p);
   static void delete_muon_pogcLcLHLT(void *p);
   static void deleteArray_muon_pogcLcLHLT(void *p);
   static void destruct_muon_pogcLcLHLT(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::muon_pog::HLT*)
   {
      ::muon_pog::HLT *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::muon_pog::HLT >(0);
      static ::ROOT::TGenericClassInfo 
         instance("muon_pog::HLT", ::muon_pog::HLT::Class_Version(), "../src/MuonPogTree.h", 313,
                  typeid(::muon_pog::HLT), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::muon_pog::HLT::Dictionary, isa_proxy, 4,
                  sizeof(::muon_pog::HLT) );
      instance.SetNew(&new_muon_pogcLcLHLT);
      instance.SetNewArray(&newArray_muon_pogcLcLHLT);
      instance.SetDelete(&delete_muon_pogcLcLHLT);
      instance.SetDeleteArray(&deleteArray_muon_pogcLcLHLT);
      instance.SetDestructor(&destruct_muon_pogcLcLHLT);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::muon_pog::HLT*)
   {
      return GenerateInitInstanceLocal((::muon_pog::HLT*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::muon_pog::HLT*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_muon_pogcLcLEventId(void *p = 0);
   static void *newArray_muon_pogcLcLEventId(Long_t size, void *p);
   static void delete_muon_pogcLcLEventId(void *p);
   static void deleteArray_muon_pogcLcLEventId(void *p);
   static void destruct_muon_pogcLcLEventId(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::muon_pog::EventId*)
   {
      ::muon_pog::EventId *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::muon_pog::EventId >(0);
      static ::ROOT::TGenericClassInfo 
         instance("muon_pog::EventId", ::muon_pog::EventId::Class_Version(), "../src/MuonPogTree.h", 338,
                  typeid(::muon_pog::EventId), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::muon_pog::EventId::Dictionary, isa_proxy, 4,
                  sizeof(::muon_pog::EventId) );
      instance.SetNew(&new_muon_pogcLcLEventId);
      instance.SetNewArray(&newArray_muon_pogcLcLEventId);
      instance.SetDelete(&delete_muon_pogcLcLEventId);
      instance.SetDeleteArray(&deleteArray_muon_pogcLcLEventId);
      instance.SetDestructor(&destruct_muon_pogcLcLEventId);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::muon_pog::EventId*)
   {
      return GenerateInitInstanceLocal((::muon_pog::EventId*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::muon_pog::EventId*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_muon_pogcLcLEvent(void *p = 0);
   static void *newArray_muon_pogcLcLEvent(Long_t size, void *p);
   static void delete_muon_pogcLcLEvent(void *p);
   static void deleteArray_muon_pogcLcLEvent(void *p);
   static void destruct_muon_pogcLcLEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::muon_pog::Event*)
   {
      ::muon_pog::Event *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::muon_pog::Event >(0);
      static ::ROOT::TGenericClassInfo 
         instance("muon_pog::Event", ::muon_pog::Event::Class_Version(), "../src/MuonPogTree.h", 351,
                  typeid(::muon_pog::Event), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::muon_pog::Event::Dictionary, isa_proxy, 4,
                  sizeof(::muon_pog::Event) );
      instance.SetNew(&new_muon_pogcLcLEvent);
      instance.SetNewArray(&newArray_muon_pogcLcLEvent);
      instance.SetDelete(&delete_muon_pogcLcLEvent);
      instance.SetDeleteArray(&deleteArray_muon_pogcLcLEvent);
      instance.SetDestructor(&destruct_muon_pogcLcLEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::muon_pog::Event*)
   {
      return GenerateInitInstanceLocal((::muon_pog::Event*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::muon_pog::Event*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace muon_pog {
//______________________________________________________________________________
atomic_TClass_ptr GenInfo::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *GenInfo::Class_Name()
{
   return "muon_pog::GenInfo";
}

//______________________________________________________________________________
const char *GenInfo::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::GenInfo*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int GenInfo::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::GenInfo*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *GenInfo::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::GenInfo*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *GenInfo::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::GenInfo*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace muon_pog
namespace muon_pog {
//______________________________________________________________________________
atomic_TClass_ptr GenParticle::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *GenParticle::Class_Name()
{
   return "muon_pog::GenParticle";
}

//______________________________________________________________________________
const char *GenParticle::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::GenParticle*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int GenParticle::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::GenParticle*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *GenParticle::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::GenParticle*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *GenParticle::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::GenParticle*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace muon_pog
namespace muon_pog {
//______________________________________________________________________________
atomic_TClass_ptr METs::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *METs::Class_Name()
{
   return "muon_pog::METs";
}

//______________________________________________________________________________
const char *METs::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::METs*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int METs::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::METs*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *METs::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::METs*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *METs::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::METs*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace muon_pog
namespace muon_pog {
//______________________________________________________________________________
atomic_TClass_ptr ChambMatch::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *ChambMatch::Class_Name()
{
   return "muon_pog::ChambMatch";
}

//______________________________________________________________________________
const char *ChambMatch::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::ChambMatch*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int ChambMatch::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::ChambMatch*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ChambMatch::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::ChambMatch*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ChambMatch::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::ChambMatch*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace muon_pog
namespace muon_pog {
//______________________________________________________________________________
atomic_TClass_ptr HitInfo::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *HitInfo::Class_Name()
{
   return "muon_pog::HitInfo";
}

//______________________________________________________________________________
const char *HitInfo::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::HitInfo*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int HitInfo::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::HitInfo*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *HitInfo::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::HitInfo*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *HitInfo::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::HitInfo*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace muon_pog
namespace muon_pog {
//______________________________________________________________________________
atomic_TClass_ptr MuonFit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MuonFit::Class_Name()
{
   return "muon_pog::MuonFit";
}

//______________________________________________________________________________
const char *MuonFit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::MuonFit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MuonFit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::MuonFit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MuonFit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::MuonFit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MuonFit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::MuonFit*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace muon_pog
namespace muon_pog {
//______________________________________________________________________________
atomic_TClass_ptr Muon::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Muon::Class_Name()
{
   return "muon_pog::Muon";
}

//______________________________________________________________________________
const char *Muon::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::Muon*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Muon::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::Muon*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Muon::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::Muon*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Muon::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::Muon*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace muon_pog
namespace muon_pog {
//______________________________________________________________________________
atomic_TClass_ptr HLTObject::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *HLTObject::Class_Name()
{
   return "muon_pog::HLTObject";
}

//______________________________________________________________________________
const char *HLTObject::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::HLTObject*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int HLTObject::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::HLTObject*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *HLTObject::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::HLTObject*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *HLTObject::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::HLTObject*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace muon_pog
namespace muon_pog {
//______________________________________________________________________________
atomic_TClass_ptr L1Muon::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *L1Muon::Class_Name()
{
   return "muon_pog::L1Muon";
}

//______________________________________________________________________________
const char *L1Muon::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::L1Muon*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int L1Muon::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::L1Muon*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *L1Muon::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::L1Muon*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *L1Muon::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::L1Muon*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace muon_pog
namespace muon_pog {
//______________________________________________________________________________
atomic_TClass_ptr HLT::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *HLT::Class_Name()
{
   return "muon_pog::HLT";
}

//______________________________________________________________________________
const char *HLT::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::HLT*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int HLT::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::HLT*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *HLT::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::HLT*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *HLT::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::HLT*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace muon_pog
namespace muon_pog {
//______________________________________________________________________________
atomic_TClass_ptr EventId::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *EventId::Class_Name()
{
   return "muon_pog::EventId";
}

//______________________________________________________________________________
const char *EventId::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::EventId*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int EventId::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::EventId*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *EventId::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::EventId*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *EventId::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::EventId*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace muon_pog
namespace muon_pog {
//______________________________________________________________________________
atomic_TClass_ptr Event::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Event::Class_Name()
{
   return "muon_pog::Event";
}

//______________________________________________________________________________
const char *Event::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::Event*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Event::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::Event*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Event::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::Event*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Event::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::muon_pog::Event*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace muon_pog
namespace muon_pog {
//______________________________________________________________________________
void GenInfo::Streamer(TBuffer &R__b)
{
   // Stream an object of class muon_pog::GenInfo.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(muon_pog::GenInfo::Class(),this);
   } else {
      R__b.WriteClassBuffer(muon_pog::GenInfo::Class(),this);
   }
}

} // namespace muon_pog
namespace ROOT {
   // Wrappers around operator new
   static void *new_muon_pogcLcLGenInfo(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::muon_pog::GenInfo : new ::muon_pog::GenInfo;
   }
   static void *newArray_muon_pogcLcLGenInfo(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::muon_pog::GenInfo[nElements] : new ::muon_pog::GenInfo[nElements];
   }
   // Wrapper around operator delete
   static void delete_muon_pogcLcLGenInfo(void *p) {
      delete ((::muon_pog::GenInfo*)p);
   }
   static void deleteArray_muon_pogcLcLGenInfo(void *p) {
      delete [] ((::muon_pog::GenInfo*)p);
   }
   static void destruct_muon_pogcLcLGenInfo(void *p) {
      typedef ::muon_pog::GenInfo current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::muon_pog::GenInfo

namespace muon_pog {
//______________________________________________________________________________
void GenParticle::Streamer(TBuffer &R__b)
{
   // Stream an object of class muon_pog::GenParticle.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(muon_pog::GenParticle::Class(),this);
   } else {
      R__b.WriteClassBuffer(muon_pog::GenParticle::Class(),this);
   }
}

} // namespace muon_pog
namespace ROOT {
   // Wrappers around operator new
   static void *new_muon_pogcLcLGenParticle(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::muon_pog::GenParticle : new ::muon_pog::GenParticle;
   }
   static void *newArray_muon_pogcLcLGenParticle(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::muon_pog::GenParticle[nElements] : new ::muon_pog::GenParticle[nElements];
   }
   // Wrapper around operator delete
   static void delete_muon_pogcLcLGenParticle(void *p) {
      delete ((::muon_pog::GenParticle*)p);
   }
   static void deleteArray_muon_pogcLcLGenParticle(void *p) {
      delete [] ((::muon_pog::GenParticle*)p);
   }
   static void destruct_muon_pogcLcLGenParticle(void *p) {
      typedef ::muon_pog::GenParticle current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::muon_pog::GenParticle

namespace muon_pog {
//______________________________________________________________________________
void METs::Streamer(TBuffer &R__b)
{
   // Stream an object of class muon_pog::METs.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(muon_pog::METs::Class(),this);
   } else {
      R__b.WriteClassBuffer(muon_pog::METs::Class(),this);
   }
}

} // namespace muon_pog
namespace ROOT {
   // Wrappers around operator new
   static void *new_muon_pogcLcLMETs(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::muon_pog::METs : new ::muon_pog::METs;
   }
   static void *newArray_muon_pogcLcLMETs(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::muon_pog::METs[nElements] : new ::muon_pog::METs[nElements];
   }
   // Wrapper around operator delete
   static void delete_muon_pogcLcLMETs(void *p) {
      delete ((::muon_pog::METs*)p);
   }
   static void deleteArray_muon_pogcLcLMETs(void *p) {
      delete [] ((::muon_pog::METs*)p);
   }
   static void destruct_muon_pogcLcLMETs(void *p) {
      typedef ::muon_pog::METs current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::muon_pog::METs

namespace muon_pog {
//______________________________________________________________________________
void ChambMatch::Streamer(TBuffer &R__b)
{
   // Stream an object of class muon_pog::ChambMatch.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(muon_pog::ChambMatch::Class(),this);
   } else {
      R__b.WriteClassBuffer(muon_pog::ChambMatch::Class(),this);
   }
}

} // namespace muon_pog
namespace ROOT {
   // Wrappers around operator new
   static void *new_muon_pogcLcLChambMatch(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::muon_pog::ChambMatch : new ::muon_pog::ChambMatch;
   }
   static void *newArray_muon_pogcLcLChambMatch(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::muon_pog::ChambMatch[nElements] : new ::muon_pog::ChambMatch[nElements];
   }
   // Wrapper around operator delete
   static void delete_muon_pogcLcLChambMatch(void *p) {
      delete ((::muon_pog::ChambMatch*)p);
   }
   static void deleteArray_muon_pogcLcLChambMatch(void *p) {
      delete [] ((::muon_pog::ChambMatch*)p);
   }
   static void destruct_muon_pogcLcLChambMatch(void *p) {
      typedef ::muon_pog::ChambMatch current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::muon_pog::ChambMatch

namespace muon_pog {
//______________________________________________________________________________
void HitInfo::Streamer(TBuffer &R__b)
{
   // Stream an object of class muon_pog::HitInfo.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(muon_pog::HitInfo::Class(),this);
   } else {
      R__b.WriteClassBuffer(muon_pog::HitInfo::Class(),this);
   }
}

} // namespace muon_pog
namespace ROOT {
   // Wrappers around operator new
   static void *new_muon_pogcLcLHitInfo(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::muon_pog::HitInfo : new ::muon_pog::HitInfo;
   }
   static void *newArray_muon_pogcLcLHitInfo(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::muon_pog::HitInfo[nElements] : new ::muon_pog::HitInfo[nElements];
   }
   // Wrapper around operator delete
   static void delete_muon_pogcLcLHitInfo(void *p) {
      delete ((::muon_pog::HitInfo*)p);
   }
   static void deleteArray_muon_pogcLcLHitInfo(void *p) {
      delete [] ((::muon_pog::HitInfo*)p);
   }
   static void destruct_muon_pogcLcLHitInfo(void *p) {
      typedef ::muon_pog::HitInfo current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::muon_pog::HitInfo

namespace muon_pog {
//______________________________________________________________________________
void MuonFit::Streamer(TBuffer &R__b)
{
   // Stream an object of class muon_pog::MuonFit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(muon_pog::MuonFit::Class(),this);
   } else {
      R__b.WriteClassBuffer(muon_pog::MuonFit::Class(),this);
   }
}

} // namespace muon_pog
namespace ROOT {
   // Wrappers around operator new
   static void *new_muon_pogcLcLMuonFit(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::muon_pog::MuonFit : new ::muon_pog::MuonFit;
   }
   static void *newArray_muon_pogcLcLMuonFit(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::muon_pog::MuonFit[nElements] : new ::muon_pog::MuonFit[nElements];
   }
   // Wrapper around operator delete
   static void delete_muon_pogcLcLMuonFit(void *p) {
      delete ((::muon_pog::MuonFit*)p);
   }
   static void deleteArray_muon_pogcLcLMuonFit(void *p) {
      delete [] ((::muon_pog::MuonFit*)p);
   }
   static void destruct_muon_pogcLcLMuonFit(void *p) {
      typedef ::muon_pog::MuonFit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::muon_pog::MuonFit

namespace muon_pog {
//______________________________________________________________________________
void Muon::Streamer(TBuffer &R__b)
{
   // Stream an object of class muon_pog::Muon.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(muon_pog::Muon::Class(),this);
   } else {
      R__b.WriteClassBuffer(muon_pog::Muon::Class(),this);
   }
}

} // namespace muon_pog
namespace ROOT {
   // Wrappers around operator new
   static void *new_muon_pogcLcLMuon(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::muon_pog::Muon : new ::muon_pog::Muon;
   }
   static void *newArray_muon_pogcLcLMuon(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::muon_pog::Muon[nElements] : new ::muon_pog::Muon[nElements];
   }
   // Wrapper around operator delete
   static void delete_muon_pogcLcLMuon(void *p) {
      delete ((::muon_pog::Muon*)p);
   }
   static void deleteArray_muon_pogcLcLMuon(void *p) {
      delete [] ((::muon_pog::Muon*)p);
   }
   static void destruct_muon_pogcLcLMuon(void *p) {
      typedef ::muon_pog::Muon current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::muon_pog::Muon

namespace muon_pog {
//______________________________________________________________________________
void HLTObject::Streamer(TBuffer &R__b)
{
   // Stream an object of class muon_pog::HLTObject.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(muon_pog::HLTObject::Class(),this);
   } else {
      R__b.WriteClassBuffer(muon_pog::HLTObject::Class(),this);
   }
}

} // namespace muon_pog
namespace ROOT {
   // Wrappers around operator new
   static void *new_muon_pogcLcLHLTObject(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::muon_pog::HLTObject : new ::muon_pog::HLTObject;
   }
   static void *newArray_muon_pogcLcLHLTObject(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::muon_pog::HLTObject[nElements] : new ::muon_pog::HLTObject[nElements];
   }
   // Wrapper around operator delete
   static void delete_muon_pogcLcLHLTObject(void *p) {
      delete ((::muon_pog::HLTObject*)p);
   }
   static void deleteArray_muon_pogcLcLHLTObject(void *p) {
      delete [] ((::muon_pog::HLTObject*)p);
   }
   static void destruct_muon_pogcLcLHLTObject(void *p) {
      typedef ::muon_pog::HLTObject current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::muon_pog::HLTObject

namespace muon_pog {
//______________________________________________________________________________
void L1Muon::Streamer(TBuffer &R__b)
{
   // Stream an object of class muon_pog::L1Muon.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(muon_pog::L1Muon::Class(),this);
   } else {
      R__b.WriteClassBuffer(muon_pog::L1Muon::Class(),this);
   }
}

} // namespace muon_pog
namespace ROOT {
   // Wrappers around operator new
   static void *new_muon_pogcLcLL1Muon(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::muon_pog::L1Muon : new ::muon_pog::L1Muon;
   }
   static void *newArray_muon_pogcLcLL1Muon(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::muon_pog::L1Muon[nElements] : new ::muon_pog::L1Muon[nElements];
   }
   // Wrapper around operator delete
   static void delete_muon_pogcLcLL1Muon(void *p) {
      delete ((::muon_pog::L1Muon*)p);
   }
   static void deleteArray_muon_pogcLcLL1Muon(void *p) {
      delete [] ((::muon_pog::L1Muon*)p);
   }
   static void destruct_muon_pogcLcLL1Muon(void *p) {
      typedef ::muon_pog::L1Muon current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::muon_pog::L1Muon

namespace muon_pog {
//______________________________________________________________________________
void HLT::Streamer(TBuffer &R__b)
{
   // Stream an object of class muon_pog::HLT.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(muon_pog::HLT::Class(),this);
   } else {
      R__b.WriteClassBuffer(muon_pog::HLT::Class(),this);
   }
}

} // namespace muon_pog
namespace ROOT {
   // Wrappers around operator new
   static void *new_muon_pogcLcLHLT(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::muon_pog::HLT : new ::muon_pog::HLT;
   }
   static void *newArray_muon_pogcLcLHLT(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::muon_pog::HLT[nElements] : new ::muon_pog::HLT[nElements];
   }
   // Wrapper around operator delete
   static void delete_muon_pogcLcLHLT(void *p) {
      delete ((::muon_pog::HLT*)p);
   }
   static void deleteArray_muon_pogcLcLHLT(void *p) {
      delete [] ((::muon_pog::HLT*)p);
   }
   static void destruct_muon_pogcLcLHLT(void *p) {
      typedef ::muon_pog::HLT current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::muon_pog::HLT

namespace muon_pog {
//______________________________________________________________________________
void EventId::Streamer(TBuffer &R__b)
{
   // Stream an object of class muon_pog::EventId.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(muon_pog::EventId::Class(),this);
   } else {
      R__b.WriteClassBuffer(muon_pog::EventId::Class(),this);
   }
}

} // namespace muon_pog
namespace ROOT {
   // Wrappers around operator new
   static void *new_muon_pogcLcLEventId(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::muon_pog::EventId : new ::muon_pog::EventId;
   }
   static void *newArray_muon_pogcLcLEventId(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::muon_pog::EventId[nElements] : new ::muon_pog::EventId[nElements];
   }
   // Wrapper around operator delete
   static void delete_muon_pogcLcLEventId(void *p) {
      delete ((::muon_pog::EventId*)p);
   }
   static void deleteArray_muon_pogcLcLEventId(void *p) {
      delete [] ((::muon_pog::EventId*)p);
   }
   static void destruct_muon_pogcLcLEventId(void *p) {
      typedef ::muon_pog::EventId current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::muon_pog::EventId

namespace muon_pog {
//______________________________________________________________________________
void Event::Streamer(TBuffer &R__b)
{
   // Stream an object of class muon_pog::Event.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(muon_pog::Event::Class(),this);
   } else {
      R__b.WriteClassBuffer(muon_pog::Event::Class(),this);
   }
}

} // namespace muon_pog
namespace ROOT {
   // Wrappers around operator new
   static void *new_muon_pogcLcLEvent(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::muon_pog::Event : new ::muon_pog::Event;
   }
   static void *newArray_muon_pogcLcLEvent(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::muon_pog::Event[nElements] : new ::muon_pog::Event[nElements];
   }
   // Wrapper around operator delete
   static void delete_muon_pogcLcLEvent(void *p) {
      delete ((::muon_pog::Event*)p);
   }
   static void deleteArray_muon_pogcLcLEvent(void *p) {
      delete [] ((::muon_pog::Event*)p);
   }
   static void destruct_muon_pogcLcLEvent(void *p) {
      typedef ::muon_pog::Event current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::muon_pog::Event

namespace ROOT {
   static TClass *vectorlEstringgR_Dictionary();
   static void vectorlEstringgR_TClassManip(TClass*);
   static void *new_vectorlEstringgR(void *p = 0);
   static void *newArray_vectorlEstringgR(Long_t size, void *p);
   static void delete_vectorlEstringgR(void *p);
   static void deleteArray_vectorlEstringgR(void *p);
   static void destruct_vectorlEstringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<string>*)
   {
      vector<string> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<string>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<string>", -2, "vector", 214,
                  typeid(vector<string>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEstringgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<string>) );
      instance.SetNew(&new_vectorlEstringgR);
      instance.SetNewArray(&newArray_vectorlEstringgR);
      instance.SetDelete(&delete_vectorlEstringgR);
      instance.SetDeleteArray(&deleteArray_vectorlEstringgR);
      instance.SetDestructor(&destruct_vectorlEstringgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<string> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<string>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEstringgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<string>*)0x0)->GetClass();
      vectorlEstringgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEstringgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEstringgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<string> : new vector<string>;
   }
   static void *newArray_vectorlEstringgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<string>[nElements] : new vector<string>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEstringgR(void *p) {
      delete ((vector<string>*)p);
   }
   static void deleteArray_vectorlEstringgR(void *p) {
      delete [] ((vector<string>*)p);
   }
   static void destruct_vectorlEstringgR(void *p) {
      typedef vector<string> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<string>

namespace ROOT {
   static TClass *vectorlEmuon_pogcLcLMuonFitgR_Dictionary();
   static void vectorlEmuon_pogcLcLMuonFitgR_TClassManip(TClass*);
   static void *new_vectorlEmuon_pogcLcLMuonFitgR(void *p = 0);
   static void *newArray_vectorlEmuon_pogcLcLMuonFitgR(Long_t size, void *p);
   static void delete_vectorlEmuon_pogcLcLMuonFitgR(void *p);
   static void deleteArray_vectorlEmuon_pogcLcLMuonFitgR(void *p);
   static void destruct_vectorlEmuon_pogcLcLMuonFitgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<muon_pog::MuonFit>*)
   {
      vector<muon_pog::MuonFit> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<muon_pog::MuonFit>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<muon_pog::MuonFit>", -2, "vector", 214,
                  typeid(vector<muon_pog::MuonFit>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEmuon_pogcLcLMuonFitgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<muon_pog::MuonFit>) );
      instance.SetNew(&new_vectorlEmuon_pogcLcLMuonFitgR);
      instance.SetNewArray(&newArray_vectorlEmuon_pogcLcLMuonFitgR);
      instance.SetDelete(&delete_vectorlEmuon_pogcLcLMuonFitgR);
      instance.SetDeleteArray(&deleteArray_vectorlEmuon_pogcLcLMuonFitgR);
      instance.SetDestructor(&destruct_vectorlEmuon_pogcLcLMuonFitgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<muon_pog::MuonFit> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<muon_pog::MuonFit>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEmuon_pogcLcLMuonFitgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<muon_pog::MuonFit>*)0x0)->GetClass();
      vectorlEmuon_pogcLcLMuonFitgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEmuon_pogcLcLMuonFitgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEmuon_pogcLcLMuonFitgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<muon_pog::MuonFit> : new vector<muon_pog::MuonFit>;
   }
   static void *newArray_vectorlEmuon_pogcLcLMuonFitgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<muon_pog::MuonFit>[nElements] : new vector<muon_pog::MuonFit>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEmuon_pogcLcLMuonFitgR(void *p) {
      delete ((vector<muon_pog::MuonFit>*)p);
   }
   static void deleteArray_vectorlEmuon_pogcLcLMuonFitgR(void *p) {
      delete [] ((vector<muon_pog::MuonFit>*)p);
   }
   static void destruct_vectorlEmuon_pogcLcLMuonFitgR(void *p) {
      typedef vector<muon_pog::MuonFit> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<muon_pog::MuonFit>

namespace ROOT {
   static TClass *vectorlEmuon_pogcLcLMuongR_Dictionary();
   static void vectorlEmuon_pogcLcLMuongR_TClassManip(TClass*);
   static void *new_vectorlEmuon_pogcLcLMuongR(void *p = 0);
   static void *newArray_vectorlEmuon_pogcLcLMuongR(Long_t size, void *p);
   static void delete_vectorlEmuon_pogcLcLMuongR(void *p);
   static void deleteArray_vectorlEmuon_pogcLcLMuongR(void *p);
   static void destruct_vectorlEmuon_pogcLcLMuongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<muon_pog::Muon>*)
   {
      vector<muon_pog::Muon> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<muon_pog::Muon>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<muon_pog::Muon>", -2, "vector", 214,
                  typeid(vector<muon_pog::Muon>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEmuon_pogcLcLMuongR_Dictionary, isa_proxy, 4,
                  sizeof(vector<muon_pog::Muon>) );
      instance.SetNew(&new_vectorlEmuon_pogcLcLMuongR);
      instance.SetNewArray(&newArray_vectorlEmuon_pogcLcLMuongR);
      instance.SetDelete(&delete_vectorlEmuon_pogcLcLMuongR);
      instance.SetDeleteArray(&deleteArray_vectorlEmuon_pogcLcLMuongR);
      instance.SetDestructor(&destruct_vectorlEmuon_pogcLcLMuongR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<muon_pog::Muon> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<muon_pog::Muon>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEmuon_pogcLcLMuongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<muon_pog::Muon>*)0x0)->GetClass();
      vectorlEmuon_pogcLcLMuongR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEmuon_pogcLcLMuongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEmuon_pogcLcLMuongR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<muon_pog::Muon> : new vector<muon_pog::Muon>;
   }
   static void *newArray_vectorlEmuon_pogcLcLMuongR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<muon_pog::Muon>[nElements] : new vector<muon_pog::Muon>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEmuon_pogcLcLMuongR(void *p) {
      delete ((vector<muon_pog::Muon>*)p);
   }
   static void deleteArray_vectorlEmuon_pogcLcLMuongR(void *p) {
      delete [] ((vector<muon_pog::Muon>*)p);
   }
   static void destruct_vectorlEmuon_pogcLcLMuongR(void *p) {
      typedef vector<muon_pog::Muon> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<muon_pog::Muon>

namespace ROOT {
   static TClass *vectorlEmuon_pogcLcLL1MuongR_Dictionary();
   static void vectorlEmuon_pogcLcLL1MuongR_TClassManip(TClass*);
   static void *new_vectorlEmuon_pogcLcLL1MuongR(void *p = 0);
   static void *newArray_vectorlEmuon_pogcLcLL1MuongR(Long_t size, void *p);
   static void delete_vectorlEmuon_pogcLcLL1MuongR(void *p);
   static void deleteArray_vectorlEmuon_pogcLcLL1MuongR(void *p);
   static void destruct_vectorlEmuon_pogcLcLL1MuongR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<muon_pog::L1Muon>*)
   {
      vector<muon_pog::L1Muon> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<muon_pog::L1Muon>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<muon_pog::L1Muon>", -2, "vector", 214,
                  typeid(vector<muon_pog::L1Muon>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEmuon_pogcLcLL1MuongR_Dictionary, isa_proxy, 0,
                  sizeof(vector<muon_pog::L1Muon>) );
      instance.SetNew(&new_vectorlEmuon_pogcLcLL1MuongR);
      instance.SetNewArray(&newArray_vectorlEmuon_pogcLcLL1MuongR);
      instance.SetDelete(&delete_vectorlEmuon_pogcLcLL1MuongR);
      instance.SetDeleteArray(&deleteArray_vectorlEmuon_pogcLcLL1MuongR);
      instance.SetDestructor(&destruct_vectorlEmuon_pogcLcLL1MuongR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<muon_pog::L1Muon> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<muon_pog::L1Muon>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEmuon_pogcLcLL1MuongR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<muon_pog::L1Muon>*)0x0)->GetClass();
      vectorlEmuon_pogcLcLL1MuongR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEmuon_pogcLcLL1MuongR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEmuon_pogcLcLL1MuongR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<muon_pog::L1Muon> : new vector<muon_pog::L1Muon>;
   }
   static void *newArray_vectorlEmuon_pogcLcLL1MuongR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<muon_pog::L1Muon>[nElements] : new vector<muon_pog::L1Muon>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEmuon_pogcLcLL1MuongR(void *p) {
      delete ((vector<muon_pog::L1Muon>*)p);
   }
   static void deleteArray_vectorlEmuon_pogcLcLL1MuongR(void *p) {
      delete [] ((vector<muon_pog::L1Muon>*)p);
   }
   static void destruct_vectorlEmuon_pogcLcLL1MuongR(void *p) {
      typedef vector<muon_pog::L1Muon> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<muon_pog::L1Muon>

namespace ROOT {
   static TClass *vectorlEmuon_pogcLcLHitInfogR_Dictionary();
   static void vectorlEmuon_pogcLcLHitInfogR_TClassManip(TClass*);
   static void *new_vectorlEmuon_pogcLcLHitInfogR(void *p = 0);
   static void *newArray_vectorlEmuon_pogcLcLHitInfogR(Long_t size, void *p);
   static void delete_vectorlEmuon_pogcLcLHitInfogR(void *p);
   static void deleteArray_vectorlEmuon_pogcLcLHitInfogR(void *p);
   static void destruct_vectorlEmuon_pogcLcLHitInfogR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<muon_pog::HitInfo>*)
   {
      vector<muon_pog::HitInfo> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<muon_pog::HitInfo>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<muon_pog::HitInfo>", -2, "vector", 214,
                  typeid(vector<muon_pog::HitInfo>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEmuon_pogcLcLHitInfogR_Dictionary, isa_proxy, 4,
                  sizeof(vector<muon_pog::HitInfo>) );
      instance.SetNew(&new_vectorlEmuon_pogcLcLHitInfogR);
      instance.SetNewArray(&newArray_vectorlEmuon_pogcLcLHitInfogR);
      instance.SetDelete(&delete_vectorlEmuon_pogcLcLHitInfogR);
      instance.SetDeleteArray(&deleteArray_vectorlEmuon_pogcLcLHitInfogR);
      instance.SetDestructor(&destruct_vectorlEmuon_pogcLcLHitInfogR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<muon_pog::HitInfo> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<muon_pog::HitInfo>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEmuon_pogcLcLHitInfogR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<muon_pog::HitInfo>*)0x0)->GetClass();
      vectorlEmuon_pogcLcLHitInfogR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEmuon_pogcLcLHitInfogR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEmuon_pogcLcLHitInfogR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<muon_pog::HitInfo> : new vector<muon_pog::HitInfo>;
   }
   static void *newArray_vectorlEmuon_pogcLcLHitInfogR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<muon_pog::HitInfo>[nElements] : new vector<muon_pog::HitInfo>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEmuon_pogcLcLHitInfogR(void *p) {
      delete ((vector<muon_pog::HitInfo>*)p);
   }
   static void deleteArray_vectorlEmuon_pogcLcLHitInfogR(void *p) {
      delete [] ((vector<muon_pog::HitInfo>*)p);
   }
   static void destruct_vectorlEmuon_pogcLcLHitInfogR(void *p) {
      typedef vector<muon_pog::HitInfo> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<muon_pog::HitInfo>

namespace ROOT {
   static TClass *vectorlEmuon_pogcLcLHLTObjectgR_Dictionary();
   static void vectorlEmuon_pogcLcLHLTObjectgR_TClassManip(TClass*);
   static void *new_vectorlEmuon_pogcLcLHLTObjectgR(void *p = 0);
   static void *newArray_vectorlEmuon_pogcLcLHLTObjectgR(Long_t size, void *p);
   static void delete_vectorlEmuon_pogcLcLHLTObjectgR(void *p);
   static void deleteArray_vectorlEmuon_pogcLcLHLTObjectgR(void *p);
   static void destruct_vectorlEmuon_pogcLcLHLTObjectgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<muon_pog::HLTObject>*)
   {
      vector<muon_pog::HLTObject> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<muon_pog::HLTObject>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<muon_pog::HLTObject>", -2, "vector", 214,
                  typeid(vector<muon_pog::HLTObject>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEmuon_pogcLcLHLTObjectgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<muon_pog::HLTObject>) );
      instance.SetNew(&new_vectorlEmuon_pogcLcLHLTObjectgR);
      instance.SetNewArray(&newArray_vectorlEmuon_pogcLcLHLTObjectgR);
      instance.SetDelete(&delete_vectorlEmuon_pogcLcLHLTObjectgR);
      instance.SetDeleteArray(&deleteArray_vectorlEmuon_pogcLcLHLTObjectgR);
      instance.SetDestructor(&destruct_vectorlEmuon_pogcLcLHLTObjectgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<muon_pog::HLTObject> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<muon_pog::HLTObject>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEmuon_pogcLcLHLTObjectgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<muon_pog::HLTObject>*)0x0)->GetClass();
      vectorlEmuon_pogcLcLHLTObjectgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEmuon_pogcLcLHLTObjectgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEmuon_pogcLcLHLTObjectgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<muon_pog::HLTObject> : new vector<muon_pog::HLTObject>;
   }
   static void *newArray_vectorlEmuon_pogcLcLHLTObjectgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<muon_pog::HLTObject>[nElements] : new vector<muon_pog::HLTObject>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEmuon_pogcLcLHLTObjectgR(void *p) {
      delete ((vector<muon_pog::HLTObject>*)p);
   }
   static void deleteArray_vectorlEmuon_pogcLcLHLTObjectgR(void *p) {
      delete [] ((vector<muon_pog::HLTObject>*)p);
   }
   static void destruct_vectorlEmuon_pogcLcLHLTObjectgR(void *p) {
      typedef vector<muon_pog::HLTObject> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<muon_pog::HLTObject>

namespace ROOT {
   static TClass *vectorlEmuon_pogcLcLGenParticlegR_Dictionary();
   static void vectorlEmuon_pogcLcLGenParticlegR_TClassManip(TClass*);
   static void *new_vectorlEmuon_pogcLcLGenParticlegR(void *p = 0);
   static void *newArray_vectorlEmuon_pogcLcLGenParticlegR(Long_t size, void *p);
   static void delete_vectorlEmuon_pogcLcLGenParticlegR(void *p);
   static void deleteArray_vectorlEmuon_pogcLcLGenParticlegR(void *p);
   static void destruct_vectorlEmuon_pogcLcLGenParticlegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<muon_pog::GenParticle>*)
   {
      vector<muon_pog::GenParticle> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<muon_pog::GenParticle>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<muon_pog::GenParticle>", -2, "vector", 214,
                  typeid(vector<muon_pog::GenParticle>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEmuon_pogcLcLGenParticlegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<muon_pog::GenParticle>) );
      instance.SetNew(&new_vectorlEmuon_pogcLcLGenParticlegR);
      instance.SetNewArray(&newArray_vectorlEmuon_pogcLcLGenParticlegR);
      instance.SetDelete(&delete_vectorlEmuon_pogcLcLGenParticlegR);
      instance.SetDeleteArray(&deleteArray_vectorlEmuon_pogcLcLGenParticlegR);
      instance.SetDestructor(&destruct_vectorlEmuon_pogcLcLGenParticlegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<muon_pog::GenParticle> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<muon_pog::GenParticle>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEmuon_pogcLcLGenParticlegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<muon_pog::GenParticle>*)0x0)->GetClass();
      vectorlEmuon_pogcLcLGenParticlegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEmuon_pogcLcLGenParticlegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEmuon_pogcLcLGenParticlegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<muon_pog::GenParticle> : new vector<muon_pog::GenParticle>;
   }
   static void *newArray_vectorlEmuon_pogcLcLGenParticlegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<muon_pog::GenParticle>[nElements] : new vector<muon_pog::GenParticle>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEmuon_pogcLcLGenParticlegR(void *p) {
      delete ((vector<muon_pog::GenParticle>*)p);
   }
   static void deleteArray_vectorlEmuon_pogcLcLGenParticlegR(void *p) {
      delete [] ((vector<muon_pog::GenParticle>*)p);
   }
   static void destruct_vectorlEmuon_pogcLcLGenParticlegR(void *p) {
      typedef vector<muon_pog::GenParticle> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<muon_pog::GenParticle>

namespace ROOT {
   static TClass *vectorlEmuon_pogcLcLGenInfogR_Dictionary();
   static void vectorlEmuon_pogcLcLGenInfogR_TClassManip(TClass*);
   static void *new_vectorlEmuon_pogcLcLGenInfogR(void *p = 0);
   static void *newArray_vectorlEmuon_pogcLcLGenInfogR(Long_t size, void *p);
   static void delete_vectorlEmuon_pogcLcLGenInfogR(void *p);
   static void deleteArray_vectorlEmuon_pogcLcLGenInfogR(void *p);
   static void destruct_vectorlEmuon_pogcLcLGenInfogR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<muon_pog::GenInfo>*)
   {
      vector<muon_pog::GenInfo> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<muon_pog::GenInfo>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<muon_pog::GenInfo>", -2, "vector", 214,
                  typeid(vector<muon_pog::GenInfo>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEmuon_pogcLcLGenInfogR_Dictionary, isa_proxy, 4,
                  sizeof(vector<muon_pog::GenInfo>) );
      instance.SetNew(&new_vectorlEmuon_pogcLcLGenInfogR);
      instance.SetNewArray(&newArray_vectorlEmuon_pogcLcLGenInfogR);
      instance.SetDelete(&delete_vectorlEmuon_pogcLcLGenInfogR);
      instance.SetDeleteArray(&deleteArray_vectorlEmuon_pogcLcLGenInfogR);
      instance.SetDestructor(&destruct_vectorlEmuon_pogcLcLGenInfogR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<muon_pog::GenInfo> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<muon_pog::GenInfo>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEmuon_pogcLcLGenInfogR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<muon_pog::GenInfo>*)0x0)->GetClass();
      vectorlEmuon_pogcLcLGenInfogR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEmuon_pogcLcLGenInfogR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEmuon_pogcLcLGenInfogR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<muon_pog::GenInfo> : new vector<muon_pog::GenInfo>;
   }
   static void *newArray_vectorlEmuon_pogcLcLGenInfogR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<muon_pog::GenInfo>[nElements] : new vector<muon_pog::GenInfo>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEmuon_pogcLcLGenInfogR(void *p) {
      delete ((vector<muon_pog::GenInfo>*)p);
   }
   static void deleteArray_vectorlEmuon_pogcLcLGenInfogR(void *p) {
      delete [] ((vector<muon_pog::GenInfo>*)p);
   }
   static void destruct_vectorlEmuon_pogcLcLGenInfogR(void *p) {
      typedef vector<muon_pog::GenInfo> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<muon_pog::GenInfo>

namespace ROOT {
   static TClass *vectorlEmuon_pogcLcLChambMatchgR_Dictionary();
   static void vectorlEmuon_pogcLcLChambMatchgR_TClassManip(TClass*);
   static void *new_vectorlEmuon_pogcLcLChambMatchgR(void *p = 0);
   static void *newArray_vectorlEmuon_pogcLcLChambMatchgR(Long_t size, void *p);
   static void delete_vectorlEmuon_pogcLcLChambMatchgR(void *p);
   static void deleteArray_vectorlEmuon_pogcLcLChambMatchgR(void *p);
   static void destruct_vectorlEmuon_pogcLcLChambMatchgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<muon_pog::ChambMatch>*)
   {
      vector<muon_pog::ChambMatch> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<muon_pog::ChambMatch>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<muon_pog::ChambMatch>", -2, "vector", 214,
                  typeid(vector<muon_pog::ChambMatch>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEmuon_pogcLcLChambMatchgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<muon_pog::ChambMatch>) );
      instance.SetNew(&new_vectorlEmuon_pogcLcLChambMatchgR);
      instance.SetNewArray(&newArray_vectorlEmuon_pogcLcLChambMatchgR);
      instance.SetDelete(&delete_vectorlEmuon_pogcLcLChambMatchgR);
      instance.SetDeleteArray(&deleteArray_vectorlEmuon_pogcLcLChambMatchgR);
      instance.SetDestructor(&destruct_vectorlEmuon_pogcLcLChambMatchgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<muon_pog::ChambMatch> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<muon_pog::ChambMatch>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEmuon_pogcLcLChambMatchgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<muon_pog::ChambMatch>*)0x0)->GetClass();
      vectorlEmuon_pogcLcLChambMatchgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEmuon_pogcLcLChambMatchgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEmuon_pogcLcLChambMatchgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<muon_pog::ChambMatch> : new vector<muon_pog::ChambMatch>;
   }
   static void *newArray_vectorlEmuon_pogcLcLChambMatchgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<muon_pog::ChambMatch>[nElements] : new vector<muon_pog::ChambMatch>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEmuon_pogcLcLChambMatchgR(void *p) {
      delete ((vector<muon_pog::ChambMatch>*)p);
   }
   static void deleteArray_vectorlEmuon_pogcLcLChambMatchgR(void *p) {
      delete [] ((vector<muon_pog::ChambMatch>*)p);
   }
   static void destruct_vectorlEmuon_pogcLcLChambMatchgR(void *p) {
      typedef vector<muon_pog::ChambMatch> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<muon_pog::ChambMatch>

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = 0);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 214,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<int>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)0x0)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace ROOT {
   static TClass *vectorlEboolgR_Dictionary();
   static void vectorlEboolgR_TClassManip(TClass*);
   static void *new_vectorlEboolgR(void *p = 0);
   static void *newArray_vectorlEboolgR(Long_t size, void *p);
   static void delete_vectorlEboolgR(void *p);
   static void deleteArray_vectorlEboolgR(void *p);
   static void destruct_vectorlEboolgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<bool>*)
   {
      vector<bool> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<bool>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<bool>", -2, "vector", 541,
                  typeid(vector<bool>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEboolgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<bool>) );
      instance.SetNew(&new_vectorlEboolgR);
      instance.SetNewArray(&newArray_vectorlEboolgR);
      instance.SetDelete(&delete_vectorlEboolgR);
      instance.SetDeleteArray(&deleteArray_vectorlEboolgR);
      instance.SetDestructor(&destruct_vectorlEboolgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<bool> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<bool>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEboolgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<bool>*)0x0)->GetClass();
      vectorlEboolgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEboolgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEboolgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<bool> : new vector<bool>;
   }
   static void *newArray_vectorlEboolgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<bool>[nElements] : new vector<bool>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEboolgR(void *p) {
      delete ((vector<bool>*)p);
   }
   static void deleteArray_vectorlEboolgR(void *p) {
      delete [] ((vector<bool>*)p);
   }
   static void destruct_vectorlEboolgR(void *p) {
      typedef vector<bool> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<bool>

namespace {
  void TriggerDictionaryInitialization_MuonPogTreeDict_Impl() {
    static const char* headers[] = {
"../src/MuonPogTree.h",
0
    };
    static const char* includePaths[] = {
"/cvmfs/cms.cern.ch/slc6_amd64_gcc530/lcg/root/6.08.07/include",
"/afs/cern.ch/work/c/calderon/private/CMSSW_9_2_1/src/MuonPOG/Tools/efficiency_comparison/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "MuonPogTreeDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace muon_pog{class __attribute__((annotate("$clingAutoload$MuonPogTree.h")))  HitInfo;}
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
namespace muon_pog{class __attribute__((annotate("$clingAutoload$MuonPogTree.h")))  ChambMatch;}
namespace muon_pog{class __attribute__((annotate("$clingAutoload$MuonPogTree.h")))  HLTObject;}
namespace muon_pog{class __attribute__((annotate("$clingAutoload$MuonPogTree.h")))  GenInfo;}
namespace muon_pog{class __attribute__((annotate("$clingAutoload$MuonPogTree.h")))  GenParticle;}
namespace muon_pog{class __attribute__((annotate("$clingAutoload$MuonPogTree.h")))  Muon;}
namespace muon_pog{class __attribute__((annotate("$clingAutoload$MuonPogTree.h")))  METs;}
namespace muon_pog{class __attribute__((annotate("$clingAutoload$MuonPogTree.h")))  MuonFit;}
namespace muon_pog{class __attribute__((annotate("$clingAutoload$MuonPogTree.h")))  L1Muon;}
namespace muon_pog{class __attribute__((annotate("$clingAutoload$MuonPogTree.h")))  HLT;}
namespace muon_pog{class __attribute__((annotate("$clingAutoload$MuonPogTree.h")))  EventId;}
namespace muon_pog{class __attribute__((annotate("$clingAutoload$MuonPogTree.h")))  Event;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "MuonPogTreeDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "../src/MuonPogTree.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"muon_pog::ChambMatch", payloadCode, "@",
"muon_pog::Event", payloadCode, "@",
"muon_pog::EventId", payloadCode, "@",
"muon_pog::GenInfo", payloadCode, "@",
"muon_pog::GenParticle", payloadCode, "@",
"muon_pog::HLT", payloadCode, "@",
"muon_pog::HLTObject", payloadCode, "@",
"muon_pog::HitInfo", payloadCode, "@",
"muon_pog::L1Muon", payloadCode, "@",
"muon_pog::METs", payloadCode, "@",
"muon_pog::Muon", payloadCode, "@",
"muon_pog::MuonFit", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("MuonPogTreeDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_MuonPogTreeDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_MuonPogTreeDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_MuonPogTreeDict() {
  TriggerDictionaryInitialization_MuonPogTreeDict_Impl();
}
